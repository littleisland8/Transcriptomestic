#Download and process GTF file in order to add chr prefix

wget https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens -o resources/Homo_sapiens.GRCh38.111.gtf.gz

bgzip -d resources/Homo_sapiens.GRCh38.111.gtf.gz

sed -i -e '/#/! s/^/chr/' resources/Homo_sapiens.GRCh38.111.gtf