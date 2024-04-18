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
    threads: 2
    wrapper:
        "v3.8.0/bio/hisat2/align"