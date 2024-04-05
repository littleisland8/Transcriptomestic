rule DESeqFromHTSeqCount:
	input:
		table="resources/example_sampleTable.tsv",
		counts=expand(f"count/htseq/{{sample}}.{{aligner}}.gene.count", sample=config["samples"].values(), aligner=config["aligners"].values())
	output:
		directory("results/STAR_HTseq")
	threads: 1
	conda:
		"../envs/DESeqFromHTSeq.yaml"
	message:
		"DESeq2 analysis of STAR HTSeq count pipeline"
	log:
		"logs/DESeqFromHTSeq.log"
	params:
		script="workflow/scripts/DESeqFromHTSeq.R",
		counts=config["DESeqFromHTSeqCount"]["count"],
		alpha=config["DESeqFromHTSeqCount"]["alpha"], 
		height=config["DESeqFromHTSeqCount"]["height"],
		width=config["DESeqFromHTSeqCount"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"
