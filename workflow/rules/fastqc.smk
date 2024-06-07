rule fastqcUntrimmed:
	input:
		config["datadir"] + "/" + "{sample}_{strand}.fq.gz"
	output:
		html="qc/{sample}_{strand}_fastqc.html",
		zip="qc/{sample}_{strand}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcUntrimmed.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"

rule fastqcTrimmed:
	input:
		"data/{sample}_{strand}.tr.fq.gz"
	output:
		html="qc/{sample}_{strand}_tr_fastqc.html",
		zip="qc/{sample}_{strand}_tr_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
	params:
		extra = "--quiet"
	log:
		"logs/{sample}_{strand}.fastqcTrimmed.log"
	threads: 1
	resources:
		mem_mb = 1024
	wrapper:
		"v3.3.3/bio/fastqc"
