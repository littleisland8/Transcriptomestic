rule featureCounts:
	input:
		bam="alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
		bai="alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai"
	output:
		"count/featureCounts/{sample}.STAR.gene.count.featureCounts.tsv"
	threads: 10
	conda:
		"../envs/featureCounts.yaml"
	message:
		"featureCounts count {wildcards.sample}"
	log:
		"logs/{sample}.featureCountsGene.log"
	params:
		stranded=config["strandedfeatureCounts"],
		fasta=config["fasta"],
		gtf=config["gtf"],
		extra=config["featureCounts_extra"]
	shell:
		"featureCounts {params.extra} -J -Q 30 -g gene_id -s {params.stranded} -G {params.fasta} -T {threads} -a {params.gtf} {input.bam} -o {output} 2>{log}"
