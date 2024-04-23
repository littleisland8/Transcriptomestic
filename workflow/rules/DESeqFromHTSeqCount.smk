rule DESeqFromHTSeqCountSTAR:
	input:
		table="resources/example_sampleTable.STAR.HTSeq.tsv",
		counts=expand(f"count/htseq/STAR/{{sample}}.STAR.gene.count", sample=config["samples"].values())
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
		script="workflow/scripts/DESeqFromHTSeq.STAR.R",
		counts=config["DESeqFromHTSeqCount"]["count"],
		alpha=config["DESeqFromHTSeqCount"]["alpha"], 
		height=config["DESeqFromHTSeqCount"]["height"],
		width=config["DESeqFromHTSeqCount"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"

rule DESeqFromHTSeqCountHisat2:
	input:
		table="resources/example_sampleTable.hisat2.HTSeq.tsv",
		counts=expand(f"count/htseq/hisat2/{{sample}}.hisat2.gene.count", sample=config["samples"].values())
	output:
		directory("results/Hisat2_HTseq")
	threads: 1
	conda:
		"../envs/DESeqFromHTSeq.yaml"
	message:
		"DESeq2 analysis of Hisat2 HTSeq count pipeline"
	log:
		"logs/DESeqFromHTSeq.log"
	params:
		script="workflow/scripts/DESeqFromHTSeq.hisat2.R",
		counts=config["DESeqFromHTSeqCount"]["count"],
		alpha=config["DESeqFromHTSeqCount"]["alpha"], 
		height=config["DESeqFromHTSeqCount"]["height"],
		width=config["DESeqFromHTSeqCount"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"
