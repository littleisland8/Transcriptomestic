rule DESeqFromFeatureCounts:
	input:
		table="resources/example_sampleTable.tsv",
		counts=expand(f"count/featureCounts/{{sample}}.{{aligner}}.STAR.gene.count.featureCounts.tsv", sample=config["samples"].values(), aligner=config["aligners"].values())
	output:
		directory("results/STAR_featureCounts/")
	threads: 1
	conda:
		"../envs/DESeqFromFeatureCounts.yaml"
	message:
		"DESeq2 analysis of STAR featureCounts count pipeline"
	log:
		"logs/DESeqFromFeatureCounts.log"
	params:
		script="workflow/scripts/DESeqFromFeatureCounts.R",
		counts=config["DESeqFromFeatureCounts"]["count"],
		alpha=config["DESeqFromFeatureCounts"]["alpha"], 
		height=config["DESeqFromFeatureCounts"]["height"],
		width=config["DESeqFromFeatureCounts"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"
