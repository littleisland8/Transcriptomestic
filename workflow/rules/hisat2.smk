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
        reads=["data/{sample}.R1.tr.fastq.gz", "data/{sample}.R2.tr.fastq.gz"],
        idx="resources" + "/" + config["genome"].split("resources/")[1].split(".fa")[0] + "_hisat2"
    output:
        "alignments/{sample}.hisat2.bam",
    log:
        "logs/{sample}.hisat2Align.log",
    params:
        extra="",
    threads: 5
    wrapper:
        "v3.8.0/bio/hisat2/align"

rule SamtoolsIndexHisat2:
    input:
        "alignments/{sample}.hisat2.bam"
    output:
        "alignments/{sample}.hisat2.bam.bai",
    log:
        "logs/{sample}.SamtoolsIndexHisat2.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"