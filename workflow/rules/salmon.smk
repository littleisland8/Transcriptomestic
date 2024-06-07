rule SalmonIndex:
	input:
		fasta=config["gentrome"]
	output:
		directory("resources/salmon_idx")
	threads: 1
	log:
		"logs/SalmonIndex.log"
	conda:
		"../envs/salmon.yaml"
	params:
		decoy="resources/decoys.txt"
	shell:
		"salmon index -t {input.fasta} -d {params.decoy} -p {threads} -i {output} 2>{log}"

rule SalmonQuant:
	input:
		idx=directory("resources/salmon_idx"),
		R1="data/{sample}_1.tr.fq.gz",
		R2="data/{sample}_2.tr.fq.gz"
	output:
		"count/{sample}_salmon/quant.sf"
	threads: 1
	log:
		"logs/{sample}.SalmonQuant.log"
	conda:
		"../envs/salmon.yaml"
	params:
		libType=config["libType"]
	shell:
		"salmon quant -p {threads} -i {input.idx} -l {params.libType} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o count/{wildcards.sample}_salmon --validateMappings --gcBias 2>{log}"


rule DESeqFromSalmon:
	input:
		table="resources/example_sampleTable.salmon.tsv",
		counts=expand(f"count/{{sample}}_salmon/quant.sf", sample=config["samples"].values())
	output:
		directory("results/Salmon/")
	threads: 1
	conda:
		"../envs/DESeqFromSalmon.yaml"
	message:
		"DESeq2 analysis of Salmon pipeline"
	log:
		"logs/DESeqFromSalmon.log"
	params:
		script="workflow/scripts/DESeqFromSalmon.R",
		counts=config["DESeqFromSalmon"]["count"],
		alpha=config["DESeqFromSalmon"]["alpha"], 
		height=config["DESeqFromSalmon"]["height"],
		width=config["DESeqFromSalmon"]["width"]
	shell:
		"Rscript {params.script} -t {input.table} -c {params.counts} -a {params.alpha} -H {params.height} -w {params.width} -o {output} 2>{log}"