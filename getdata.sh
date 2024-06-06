#get fasta 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O resources/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd resources

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

cd -

#Download and process GTF file and add chr prefix

wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz -P resources

bgzip -d resources/Homo_sapiens.GRCh38.111.gtf.gz

sed -i -e '/#/! s/^/chr/' resources/Homo_sapiens.GRCh38.111.gtf

#Make gentrome 

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P resources/

cd resources 

grep "^>" GRCh38_full_analysis_set_plus_decoy_hla.fa | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

bgzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz

cat Homo_sapiens.GRCh38.cdna.all.fa GRCh38_full_analysis_set_plus_decoy_hla.fa > gentrome.fa && bgzip gentrome.fa

cd -