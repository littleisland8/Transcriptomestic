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
        temp("alignments/{sample}.STAR.Aligned.out.sam")
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
        "STAR --runThreadN {threads} --readFilesCommand {params.readFilesCommand} --outFileNamePrefix alignments/{wildcards.sample}.STAR. --genomeDir {input.idx} --readFilesIn {input.R1} {input.R2} 2>{log}"

rule SamtoolsSort:
    input:
        "alignments/{sample}.STAR.Aligned.out.sam"
    output:
        "alignments/{sample}.STAR.srt.bam",
    log:
       "logs/{sample}.SamtoolsSort.log"
    params:
        extra="-m 4G",
    threads: 8
    wrapper:
        "v3.3.6/bio/samtools/sort"

rule SamtoolsIndex:
    input:
        "alignments/{sample}.STAR.srt.bam",
    output:
        "alignments/{sample}.STAR.srt.bam.bai",
    log:
        "logs/{sample}.SamtoolsIndex.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"