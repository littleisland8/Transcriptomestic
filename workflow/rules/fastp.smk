rule fastpTumor: ## aggiungere flag per i threads
	input:
		sample=["data/{sample}.R1.fastq.gz", "data/{sample}.R2.fastq.gz"]
	output:
		trimmed=["data/{sample}.R1.tr.fastq.gz", "data/{sample}.R2.tr.fastq.gz"],
		json="data/{sample}.json",
		failed="data/{sample}.failedreads.txt",
		html="data/{sample}.html",
		unpaired1="data/{sample}.u1.fq.gz",
		unpaired2="data/{sample}.u2.fq.gz"
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