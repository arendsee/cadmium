#!/bin/bash
set -u

base=$PWD
db=$PWD/brassicaceae

# Arabidopsis_thaliana
at=$db/Arabidopsis_thaliana
at_fna=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.fna.gz
at_gff=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_genomic.gff.gz
at_pep=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_protein.faa.gz
# at_trs=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001735.3_TAIR10/GCF_000001735.3_TAIR10_rna.fna.gz
wget -O ${at}.fna ${at_fna}.gz
wget -O ${at}.gff ${at_gff}.gz
wget -O ${at}.faa ${at_pep}.gz
# wget -O ${at}.transcript.fna ${at_trs}.gz

# Arabidopsis_lyrata
al=$db/Arabidopsis_lyrata
al_fna=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_genomic.fna.gz
al_gff=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_genomic.gff.gz
al_pep=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_protein.faa.gz
# al_trs=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000004255.1_v.1.0/GCF_000004255.1_v.1.0_rna.fna.gz
wget -O ${al}.fna ${al_fna}.gz
wget -O ${al}.gff ${al_gff}.gz
wget -O ${al}.faa ${al_pep}.gz
# wget -O ${al}.transcript.fna ${al_trs}.gz

# Capsella_rubella
cr=$db/Capsella_rubella
cr_fna=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_genomic.fna.gz
cr_gff=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_genomic.gff.gz
cr_pep=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_protein.faa.gz
# cr_trs=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_rna.fna.gz
wget -O ${cr}.fna ${cr_fna}.gz
wget -O ${cr}.gff ${cr_gff}.gz
wget -O ${cr}.faa ${cr_pep}.gz
# wget -O ${cr}.transcript.fna ${cr_trs}.gz

# Brassica_rapa
br=$db/Brassica_rapa
br_fna=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000309985.1_Brapa_1.0/GCF_000309985.1_Brapa_1.0_genomic.fna.gz
br_gff=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000309985.1_Brapa_1.0/GCF_000309985.1_Brapa_1.0_genomic.gff.gz
br_pep=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000309985.1_Brapa_1.0/GCF_000309985.1_Brapa_1.0_protein.faa.gz
# br_trs=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000309985.1_Brapa_1.0/GCF_000309985.1_Brapa_1.0_rna.fna.gz
wget -O ${br}.fna ${br_fna}.gz
wget -O ${br}.gff ${br_gff}.gz
wget -O ${br}.faa ${br_pep}.gz
# wget -O ${br}.transcript.fna ${br_trs}.gz

# Brassicaea family tree
echo '(Brassica_rapa,(Capsella_rubella,(Arabidopsis_lyrata,Arabidopsis_thaliana)))' > $db/brassicaceae.tree
