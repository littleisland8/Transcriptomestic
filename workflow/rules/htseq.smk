rule HtseqcountGeneFromSTAR:
	input:
		bam="alignments/{sample}.STAR.srt.bam",
		bai="alignments/{sample}.STAR.srt.bam.bai"
	output:
		"count/{sample}.STAR.gene.count"
	threads: 1
	conda:
		"../envs/htseq.yaml"
	message:
		"htseq count {wildcards.sample}"
	log:
		"logs/{sample}.HtseqcountGene.log"
	params:
		stranded=config["stranded"],
		mode=config["htseqmode"],
		gtf=config["gtf"]
	shell:
		"htseq-count -f bam --max-reads-in-buffer 80000000 -m {params.mode} -t gene --idattr gene_id -r pos -s {params.stranded} {input.bam} {params.gtf} > {output} 2>{log}"








		