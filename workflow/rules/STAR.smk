rule STARIndex:
    input:
        fasta=config["genome"],
    output:
        log="resources/GRCh38_full_analysis_set_plus_decoy_hla/Log.out",
        idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla")
    message:
        "STAR index"
    threads: 5
    conda:
        "../envs/STAR.yaml"
    params:
        sjdbOverhang=str(config["sjdbOverhang"]),
        gtf=config["gtf"]
    log:
        "logs/star_index_GRCh38_full_analysis_set_plus_decoy_hla.fa.log"
    shell:
        "STAR --runMode genomeGenerate --genomeDir {output.idx} --genomeFastaFiles {input} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.sjdbOverhang} 2>{log}"    


rule STARAlign:
    input:
        log="resources/GRCh38_full_analysis_set_plus_decoy_hla/Log.out",
        idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla"),
        R1="data/{sample}.R1.tr.fastq.gz",
        R2="data/{sample}.R2.tr.fastq.gz"
    output:
        "alignments/{sample}Aligned.out.sam"
    message:
        "STAR Align"
    threads: 5
    conda:
        "../envs/STAR.yaml"
    params:
        readFilesCommand="zcat",
        gtf=config["gtf"]
    log:
        "logs/{sample}.STARAlign.log"
    shell:
        "STAR --runThreadN {threads} --readFilesCommand {params.readFilesCommand} --outFileNamePrefix alignments/{wildcards.sample} --genomeDir {input.idx} --readFilesIn {input.R1} {input.R2} 2>{log}"



#rule star_index:
#    input:
#        fasta=config["genome"],
#    output:
#        directory("resources/GRCh38_full_analysis_set_plus_decoy_hla")
#    message:
#        "Testing STAR index"
#    threads: 10
#    params:
#        extra="--sjdbOverhang " + str(config["sjdbOverhang"]) + " --sjdbGTFfile " + config["gtf"]
#    log:
#        "logs/star_index_GRCh38_full_analysis_set_plus_decoy_hla.fa.log"
#    wrapper:
#        "v3.3.6/bio/star/index"

#rule StarAlign:
#    input:
        # use a list for multiple fastq files for one sample
        # usually technical replicates across lanes/flowcells
#        fq1=["data/{sample}.R1.tr.fastq.gz", "data/{sample}.R2.tr.fastq.gz"],
        # paired end reads needs to be ordered so each item in the two lists match
        #fq2=["reads/{sample}_R2.1.fastq", "reads/{sample}_R2.2.fastq"],  #optional
        # path to STAR reference genome index
#        idx="resources/GRCh38_full_analysis_set_plus_decoy_hla",
#    output:
        # see STAR manual for additional output files
#        aln=temp("alignments/{sample}.Aligned.out.sam"),
#       log="logs/{sample}.Log.out",
#        sj="alignments/{sample}.SJ.out.tab",
#        finalLog="alignments/{sample}Log.final.out"
        #unmapped=["alignments/{sample}/unmapped.1.fastq.gz","alignments/{sample}/unmapped.2.fastq.gz"],
#    log:
#        "logs/{sample}.log",
#    params:
        # optional parameters
#        extra="",
#    threads: 8
#    wrapper:
#        "v3.3.6/bio/star/align"

rule SamtoolsSort:
    input:
        "alignments/{sample}Aligned.out.sam"
    output:
        "alignments/{sample}.Aligned.out.srt.bam",
    log:
       "logs/{sample}.SamtoolsSort.log"
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v3.3.6/bio/samtools/sort"

rule SamtoolsIndex:
    input:
        "alignments/{sample}.Aligned.out.srt.bam",
    output:
        "alignments/{sample}.Aligned.out.srt.bam.bai",
    log:
        "logs/{sample}.SamtoolsIndex.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"