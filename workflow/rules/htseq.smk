rule HtseqcountGeneFromSTAR:
	input:
		bam="alignments/{sample}.{aligner}.Aligned.sortedByCoord.out.bam",
		bai="alignments/{sample}.{aligner}.Aligned.sortedByCoord.out.bam.bai"
	output:
		"count/htseq/{sample}.{aligner}.gene.count"
	threads: 1
	conda:
		"../envs/htseq.yaml"
	message:
		"htseq count {wildcards.sample}"
	log:
		"logs/{sample}.{aligner}.HtseqcountGene.log"
	params:
		stranded=config["stranded"],
		mode=config["htseqmode"],
		gtf=config["gtf"],
		extra=config["HTSeq_extra"]
	shell:
		"htseq-count -f bam --max-reads-in-buffer 80000000 -m {params.mode} -a {params.extra} -t gene --idattr gene_id -r pos -s {params.stranded} {input.bam} {params.gtf} > {output} 2>{log}"








		