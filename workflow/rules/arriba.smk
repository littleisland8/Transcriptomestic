rule arriba:
	input:
		log="resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR/Log.out",
		idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR"),
		R1=config["pipedir"] + "/" + "data/{sample}_1.tr.fq.gz",
		R2=config["pipedir"] + "/" + "data/{sample}_2.tr.fq.gz"
	output:
		fusion="fusion/{sample}.fusions.tsv",
		discarded="fusion/{sample}.fusions.discarded.tsv"
	threads: 8		
	conda:
		"../envs/arriba.yaml"
	log:
		"logs/{sample}.arriba.log"
	params:
		gtf=config["gtf"],
		gff3=config["protein_domain"],
		genome=config["genome"],
		blacklist=config["blacklist_arriba"],
		known_fusion=config["known_fusion"],
		strandedness=config["strandedAribba"]
	shell:
		"STAR --runThreadN 8 \
	--genomeDir {input.idx} --genomeLoad NoSharedMemory \
	--readFilesIn {input.R1} {input.R2} --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |arriba \
	-x /dev/stdin \
	-o {output.fusion} -O {output.discarded} \
	-a {params.genome} -g {params.gtf} \
	-k {params.known_fusion} -t {params.known_fusion} -b {params.blacklist} -p {params.gff3} -s {params.strandedness} 2>{log}" # -k /path/to/known_fusions.tsv.gz -t /path/to/known_fusions.tsv.gz -b /path/to/blacklist.tsv.gz


