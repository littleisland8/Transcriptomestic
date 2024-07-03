rule GenerateBigWigSTAR:
    input:
        bam=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
        bai=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai"
    output:
        "alignments/{sample}.STAR.bw",
    log:
        "logs/{sample}.GenerateBigWigSTAR.log",
    threads: 4  
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output} --normalizeUsing RPKM --binSize 1 --numberOfProcessors {threads} -v 2>{log}"

rule GenerateBigWigHisat2:
    input:
        bam=config["pipedir"] + "/" + "alignments/{sample}.hisat2.srt.bam",
        bai=config["pipedir"] + "/" + "alignments/{sample}.hisat2.srt.bam.bai"
    output:
        "alignments/{sample}.hisat2.bw",
    log:
        "logs/{sample}.GenerateBigWigHisat2.log",
    threads: 4  
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output} --normalizeUsing RPKM --binSize 1 --numberOfProcessors {threads} -v 2>{log}"