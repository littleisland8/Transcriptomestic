rule SalmonIndex:
	input:
		fasta=config["genome"]
	output:
		directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_salmon")
	threads: 1
	log:
		"logs/SalmonIndex.log"
	conda:
		"../envs/salmon.yaml"
	shell:
		"salmon index -t {input.fasta} -i {output} 2>{log}"

rule SalmonQuant:
	input:
		idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_salmon"),
		R1="data/{sample}.R1.tr.fastq.gz",
		R2="data/{sample}.R2.tr.fastq.gz"
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
		"salmon quant -p 20 -i {input.idx} -l {params.libType} -1 <(zcat {input.R1}) -2 <(zcat {input.R2}) -o count/{wildcards.sample}_salmon --validateMappings --gcBias 2>{log}"
