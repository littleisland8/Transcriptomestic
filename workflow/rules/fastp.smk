rule fastpTumor: ## aggiungere flag per i threads
	input:
		sample=[config["datadir"] + "/" + "{sample}_R1.fastq.gz", config["datadir"] + "/" + "{sample}_R2.fastq.gz"]
	output:
		trimmed=[config["pipedir"] + "/" + "data/{sample}_R1.tr.fastq.gz", config["pipedir"] + "/" + "data/{sample}_R2.tr.fastq.gz"],
		json=config["pipedir"] + "/" + "data/{sample}.json",
		failed=config["pipedir"] + "/" + "data/{sample}.failedreads.txt",
		html=config["pipedir"] + "/" + "data/{sample}.html",
		unpaired1=config["pipedir"] + "/" + "data/{sample}.u1.fastq.gz",
		unpaired2=config["pipedir"] + "/" + "data/{sample}.u2.fastq.gz"
	threads: 1
	log:
		"logs/{sample}.fastp.log"
	message:
		"Trimming with fastp"	
	params:
		#adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		extra=" --average_qual 25 --low_complexity_filter --complexity_threshold 30 -x --detect_adapter_for_pe --overrepresentation_analysis"
	wrapper:
		"v3.3.3/bio/fastp"
