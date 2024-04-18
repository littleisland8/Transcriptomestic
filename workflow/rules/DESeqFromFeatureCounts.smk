rule DESeqFromFeatureCountsSTAR:
	input:
		table="resources/example_sampleTable.STAR.FeatureCounts.tsv",
		counts=expand(f"count/featureCounts/{{sample}}.STAR.gene.count.featureCounts.tsv", sample=config["samples"].values())
	output:
		directory("results/STAR_featureCounts/")
	threads: 1
	conda:
		"../envs/DESeqFromFeatureCounts.yaml"
	message:
		"DESeq2 analysis of STAR featureCounts count pipeline"
	log:
		"logs/DESeqFromFeatureCountsSTAR.log"
	params:
		script="workflow/scripts/DESeqFromFeatureCounts.R",
		counts=config["DESeqFromFeatureCounts"]["count"],
		alpha=config["DESeqFromFeatureCounts"]["alpha"], 
		height=config["DESeqFromFeatureCounts"]["height"],
		width=config["DESeqFromFeatureCounts"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"

rule DESeqFromFeatureCountsHisat2:
	input:
		table="resources/example_sampleTable.hisat2.FeatureCounts.tsv",
		counts=expand(f"count/featureCounts/{{sample}}.hisat2.gene.count.featureCounts.tsv", sample=config["samples"].values())
	output:
		directory("results/Hisat2_featureCounts/")
	threads: 1
	conda:
		"../envs/DESeqFromFeatureCounts.yaml"
	message:
		"DESeq2 analysis of Hisat2 featureCounts count pipeline"
	log:
		"logs/DESeqFromFeatureCountsHisat2.log"
	params:
		script="workflow/scripts/DESeqFromFeatureCounts.R",
		counts=config["DESeqFromFeatureCounts"]["count"],
		alpha=config["DESeqFromFeatureCounts"]["alpha"], 
		height=config["DESeqFromFeatureCounts"]["height"],
		width=config["DESeqFromFeatureCounts"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"