rule GenerateBigWig:
    input:
        "alignments/{sample}.STAR.srt.bam",
    output:
        "alignments/{sample}.STAR.bw",
    log:
        "logs/{sample}.GenerateBigWig.log",
    threads: 4  
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input} -o {output} --normalizeUsing RPKM --binSize 1 --numberOfProcessors {threads} -v 2>{log}"