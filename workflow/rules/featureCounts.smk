rule featureCountsSTAR:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai"
	output:
		"count/featureCounts/STAR/{sample}.STAR.gene.count.featureCounts.tsv"
	threads: 10
	conda:
		"../envs/featureCounts.yaml"
	message:
		"featureCounts count {wildcards.sample}"
	log:
		"logs/{sample}.STAR.featureCountsSTAR.log"
	params:
		stranded=config["strandedfeatureCounts"],
		fasta=config["genome"],
		gtf=config["gtf"],
		extra=config["featureCounts_extra"]
	shell:
		"featureCounts {params.extra} -J -Q 30 -g gene_id -s {params.stranded} -G {params.fasta} -T {threads} -a {params.gtf} {input.bam} -o {output} 2>{log}"

rule featureCountsHisat2:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.hisat2.srt.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.hisat2.srt.bam.bai"
	output:
		"count/featureCounts/hisat2/{sample}.hisat2.gene.count.featureCounts.tsv"
	threads: 10
	conda:
		"../envs/featureCounts.yaml"
	message:
		"featureCounts count {wildcards.sample}"
	log:
		"logs/{sample}.hisat2.featureCountsGene.log"
	params:
		stranded=config["strandedfeatureCounts"],
		fasta=config["genome"],
		gtf=config["gtf"],
		extra=config["featureCounts_extra"]
	shell:
		"featureCounts {params.extra} -J -Q 30 -g gene_id -s {params.stranded} -G {params.fasta} -T {threads} -a {params.gtf} {input.bam} -o {output} 2>{log}"
