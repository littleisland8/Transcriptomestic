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
		strandedness=config["strandedAribba"],
		outdir=config["pipedir"] + "/" + "alignments"
	shell:
		"STAR --runThreadN 8 \
	--genomeDir {input.idx} --genomeLoad NoSharedMemory \
	--readFilesIn {input.R1} {input.R2} --outFileNamePrefix {params.outdir}/{wildcards.sample}.STAR_arriba. --readFilesCommand zcat \
	--outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 |arriba \
	-x /dev/stdin \
	-o {output.fusion} -O {output.discarded} \
	-a {params.genome} -g {params.gtf} \
	-k {params.known_fusion} -t {params.known_fusion} -b {params.blacklist} -p {params.gff3} -s {params.strandedness} 2>{log}" # -k /path/to/known_fusions.tsv.gz -t /path/to/known_fusions.tsv.gz -b /path/to/blacklist.tsv.gz

#rule STARArriba:
#	input:
#		log="resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR/Log.out",
#		idx=directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_STAR"),
#		R1=config["pipedir"] + "/" + "data/{sample}_1.tr.fq.gz",
#		R2=config["pipedir"] + "/" + "data/{sample}_2.tr.fq.gz"
#	output:
#		fusion=config["pipedir"] + "/" + "fusion/{sample}.fusions.tsv",
#		discarded=config["pipedir"] + "/" + "fusion/{sample}.fusions.discarded.tsv"
#	threads: 8		
#	conda:
#		"../envs/STAR.yaml"
#	log:
#		"logs/{sample}.ArribaSTAR.log"
#	params:
#        readFilesCommand="zcat",
#        gtf=config["gtf"],
#        outdir=config["pipedir"] + "/" + "fusion"
#	shell:
#		"STAR --runThreadN 8 \
#	--genomeDir {input.idx} --genomeLoad NoSharedMemory \
#	--readFilesIn {input.R1} {input.R2} --readFilesCommand {params.readFilesCommand} \
#	--outStd BAM_Unsorted --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMcompression 0 \
#	--outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
#	--chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
#	--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 --outFileNamePrefix {params.outdir}/{wildcards.sample}.STAR"


rule ConvertToVcf:
	input:
		"fusion/{sample}.fusions.tsv"
	output:
		vcf="fusion/{sample}.fusions.srt.vcf.gz",
		tbi="fusion/{sample}.fusions.srt.vcf.gz.tbi"
	threads: 1
	log:
		"logs/{sample}.ConvertToVcf.log"
	params:
		script="workflow/scripts/convert_fusions_to_vcf.sh",
		genome=config["genome"],
		tmp="fusion/{sample}.fusions.vcf",
		tmp2="fusion/{sample}.fusions.srt.vcf.gz"
	shell:
		"bash {params.script} {params.genome} {input} {params.tmp} 2>{log} && bcftools sort -O z -o {params.tmp2} {params.tmp} 2>>{log} && tabix -p vcf {params.tmp2} 2>>{log} && rm {params.tmp}"  

rule PairedVcf:
	input:
		expand(f"fusion/{{sample}}.fusions.srt.vcf.gz", sample=config["samples"].values())
	output:
		expand(f"fusion/{{sample}}.lactate.fusions.vcf.gz", sample=config["paired"].values())
	threads: 1
	log:
		"logs/Arriba.IntersectVcf.log"
	conda:
		"../envs/bcftools.yaml"
	params:
		sample=' '.join(list(config['paired'].values()))
	shell:
		'''
		for s in {params.sample}; do ctr=$(echo ${{s}} |cut -d "." -f1) && lact=$(echo ${{s}} |cut -d "." -f2) && bcftools isec -w 2 -n~01 -O z -o fusion/${{ctr}}.${{lact}}.lactate.fusions.vcf.gz fusion/${{ctr}}.fusions.srt.vcf.gz fusion/${{lact}}.fusions.srt.vcf.gz; done  2>>{log}
		'''

rule ParseVcf:
	input:
		"fusion/{sample}.lactate.fusions.vcf.gz"
	output:
		"fusion/{sample}.lactate.fusions.tsv"
	threads:1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.ParseVcf.log"
	params:
		header="resources/arriba.h"
	shell:
		'''
		lact=$(echo {wildcards.sample} |cut -d "." -f2) && bcftools query -f "%POS\n" {input} |grep -v -w "^1" > fusion/{wildcards.sample}.pos.txt && cat {params.header} <(grep -w -f fusion/{wildcards.sample}.pos.txt fusion/${{lact}}.fusions.tsv) > {output} 2>>{log}
		'''

rule MergeVcf:
	input:
		expand(f"fusion/{{sample}}.lactate.fusions.vcf.gz", sample=config["samples"].values())
	output:
		"fusion/total.lactate.fusion.vcf.gz"
	threads:1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}.ParseVcf.log"
	params:
		header="resources/arriba.h"
	shell:
		'''
		bcftools isec -n+2 -O z -o {output} {input}
		'''