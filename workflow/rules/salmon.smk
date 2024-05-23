rule SalmonIndex:
	input:
		fasta=config["genome"]
	output:
		directory("resources/GRCh38_full_analysis_set_plus_decoy_hla_salmon")
	threads: 1
	conda:
		"../envs/salmon.yaml"
	shell:
		"salmon index -t {input.fasta} -i {output}"