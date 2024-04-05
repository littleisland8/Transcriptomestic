rule featureCounts:
	input:
		bam="alignments/{sample}.{aligner}.Aligned.sortedByCoord.out.bam",
		bai="alignments/{sample}.{aligner}.Aligned.sortedByCoord.out.bam.bai"
	output:
		"count/featureCounts/{sample}.{aligner}.gene.count.featureCounts.tsv"
	threads: 10
	conda:
		"../envs/featureCounts.yaml"
	message:
		"featureCounts count {wildcards.sample}"
	log:
		"logs/{sample}.{aligner}.featureCountsGene.log"
	params:
		stranded=config["strandedfeatureCounts"],
		fasta=config["genome"],
		gtf=config["gtf"],
		extra=config["featureCounts_extra"]
	shell:
		"featureCounts {params.extra} -J -Q 30 -g gene_id -s {params.stranded} -G {params.fasta} -T {threads} -a {params.gtf} {input.bam} -o {output} 2>{log}"
