rule GenerateBigWig:
    input:
        bam="alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
        bai="alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai"
    output:
        "alignments/{sample}.STAR.bw",
    log:
        "logs/{sample}.GenerateBigWig.log",
    threads: 4  
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output} --normalizeUsing RPKM --binSize 1 --numberOfProcessors {threads} -v 2>{log}"