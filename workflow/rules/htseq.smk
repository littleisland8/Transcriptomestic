rule HtseqcountGeneFromSTAR:
	input:
		bam=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam",
		bai=config["pipedir"] + "/" + "alignments/{sample}.STAR.Aligned.sortedByCoord.out.bam.bai"
	output:
		"count/htseq/STAR/{sample}.STAR.gene.count"
	threads: 1
	conda:
		"../envs/htseq.yaml"
	message:
		"htseq count {wildcards.sample}"
	log:
		"logs/{sample}.HtseqcountGeneFromSTAR.log"
	params:
		stranded=config["stranded"],
		mode=config["htseqmode"],
		gtf=config["gtf"],
		extra=config["HTSeq_extra"]
	shell:
		"htseq-count -f bam --max-reads-in-buffer 80000000 -m {params.mode} -a {params.extra} -t gene --idattr gene_id -r pos -s {params.stranded} {input.bam} {params.gtf} > {output} 2>{log}"

rule HtseqcountGeneFromHisat2:
	input:
		bam="alignments/{sample}.hisat2.srt.bam",
		bai="alignments/{sample}.hisat2.srt.bam.bai"
	output:
		"count/htseq/hisat2/{sample}.hisat2.gene.count"
	threads: 1
	conda:
		"../envs/htseq.yaml"
	message:
		"htseq count {wildcards.sample}"
	log:
		"logs/{sample}.HtseqcountGeneFromHisat2.log"
	params:
		stranded=config["stranded"],
		mode=config["htseqmode"],
		gtf=config["gtf"],
		extra=config["HTSeq_extra"]
	shell:
		"htseq-count -f bam --max-reads-in-buffer 80000000 -m {params.mode} -a {params.extra} -t gene --idattr gene_id -r pos -s {params.stranded} {input.bam} {params.gtf} > {output} 2>{log}"






		
