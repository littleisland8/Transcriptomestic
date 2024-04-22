rule hisat2_index:
    input:
        fasta=config["genome"]
    output:
        directory("resources" + "/" + config["genome"].split("resources/")[1].split(".fa")[0] + "_hisat2")
    params:
        prefix = "resources" + "/" + config["genome"].split("resources/")[1].split(".fa")[0] + "_hisat2"
    log:
        "logs/hisat2_index.log"
    threads: 2
    wrapper:
        "v3.8.0/bio/hisat2/index"

rule hisat2Align:
    input:
        R1="data/{sample}.R1.tr.fastq.gz", 
        R2="data/{sample}.R2.tr.fastq.gz",
        idx="resources" + "/" + config["genome"].split("resources/")[1].split(".fa")[0] + "_hisat2"
    output:
        "alignments/{sample}.hisat2.bam"
    threads: 5
    conda: 
        "../envs/hisat2.yaml"
    log:
        "logs/{sample}.hisat2Align.log"
    params:
        strandness=config["rna_strandness"]
    shell:
        "hisat2 -p {threads} -x {input.idx} --rna-strandness {params.strandness} -1 {input.R1} -2 {input.R2} | samtools view -Sbh -o {output} 2>{log}"

rule SamtoolsSortHisat2:
    input:
        "alignments/{sample}.hisat2.bam"
    output:
        "alignments/{sample}.hisat2.srt.bam"
    log:
        "{sample}.log",
    params:
        extra="-m 4G",
    threads: 5
    wrapper:
        "v3.8.0/bio/samtools/sort"

rule SamtoolsIndexHisat2:
    input:
        "alignments/{sample}.hisat2.srt.bam"
    output:
        "alignments/{sample}.hisat2.srt.bam.bai",
    log:
        "logs/{sample}.SamtoolsIndexHisat2.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"