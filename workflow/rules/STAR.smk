rule STARIndex:
    input:
        fasta=config["genome"],
    output:
        log="resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR/Log.out",
        idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR")
    message:
        "STAR index"
    threads: 5
    conda:
        "../envs/STAR.yaml"
    params:
        sjdbOverhang=str(config["sjdbOverhang"]),
        gtf=config["gtf"]
    log:
        "logs/star_index.log"
    shell:
        "STAR --runMode genomeGenerate --genomeDir {output.idx} --genomeFastaFiles {input} --sjdbGTFfile {params.gtf} --sjdbOverhang {params.sjdbOverhang} 2>{log}"    

rule STARAlign:
    input:
        log="resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR/Log.out",
        idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR"),
        R1="data/{sample}.R1.tr.fastq.gz",
        R2="data/{sample}.R2.tr.fastq.gz"
    output:
        "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam"
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
        "STAR --runThreadN {threads} --readFilesCommand {params.readFilesCommand} --outFileNamePrefix alignments/{wildcards.sample}.STAR. --outSAMunmapped Within --outFilterMismatchNmax 10 --quantMode TranscriptomeSAM GeneCounts --chimSegmentMin 15 --outSAMtype BAM SortedByCoordinate --genomeDir {input.idx} --readFilesIn {input.R1} {input.R2} 2>{log}"

rule SamtoolsIndex:
    input:
        "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
    output:
        "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/{sample}.SamtoolsIndex.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"