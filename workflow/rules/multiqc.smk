rule multiQC_step1:
	input:
		html=expand(config["pipedir"] + "/" + f"fastqc/{{sample}}_{{strand}}_fastqc.html", sample=config["samples"].values(),strand=config["strand"].values()),
		zip_=expand(config["pipedir"] + "/" + f"fastqc/{{sample}}_{{strand}}_fastqc.zip", sample=config["samples"].values(),strand=config["strand"].values()),
		bam=expand(config["pipedir"] + "/" + f"alignments/{{sample}}.STAR.Aligned.sortedByCoord.out.bam", sample=config["samples"].values()),
		counts=expand(f"count/htseq/STAR/{{sample}}.STAR.gene.count", sample=config["samples"].values())
	output:
		"qc/multiqc_report.html"
	message:
		"MultiQC Fastq and STAR"
	threads: 5
	conda:
		"../envs/multiqc.yaml"
	log:
		"logs/multiQC_step1.log"
	params:
		pipedir=config["pipedir"]
	shell:
		'''
		multiqc {params.pipedir} count/ --outdir qc 
		'''