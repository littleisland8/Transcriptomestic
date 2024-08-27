#get fasta 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -O resources/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd resources

samtools faidx GRCh38_full_analysis_set_plus_decoy_hla.fa

cd -

#Download and process GTF file and add chr prefix

wget https://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz -P resources

bgzip -d resources/Homo_sapiens.GRCh38.112.gtf.gz

sed -i -e '/#/! s/^/chr/' resources/Homo_sapiens.GRCh38.112.gtf

#Download and process gff3 file

wget https://ftp.ensembl.org/pub/release-112/gff3/homo_sapiens/Homo_sapiens.GRCh38.112.gff3.gz -P resources

bgzip -d Homo_sapiens.GRCh38.112.gff3.gz

sed -i -e '/#/! s/^/chr/' resources/Homo_sapiens.GRCh38.112.gff3

#Make gentrome 

wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -P resources/
wget https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz

cd resources 

zgrep "^>" Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

bgzip -d Homo_sapiens.GRCh38.cdna.all.fa.gz

cat  Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz > gentrome.fa.gz

cd -

wget https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz -O arriba_v2.4.0.tar.gz -P resources/

cd resources/

tar -xzvf arriba_v2.4.0.tar.gz

rm arriba_v2.4.0.tar.gz

cd -




