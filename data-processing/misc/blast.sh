#!/bin/bash

# load the BLAST+ module
ml BLAST+/2.10.1-gompi-2022a

# unzip the Arabidopsis thaliana reference genome
echo "unzipping Araport11 protein file..."
gunzip Araport11_pep_20220914_representative_gene_model.gz

# make BLAST 
echo "creating BLAST database..."
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model -dbtype prot

# blastx of til nuc to arabidopsis thaliana protein, selecting top hit
echo "running blastx: tilingii -> Arabidopsis thaliana..."
blastx -query genes_of_interest.fasta \
       -db /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model \
       -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 \
       -out til_to_AT.txt

# tblastn of arabidopsis thaliana protein to tilingii nuc, selecting top hit
echo "running tblastn: Arabidopsis thaliana -> tilingii..."
tblastn -query /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model \
        -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 \
        -out AT_to_til.txt

# Determine reciprocal hits
echo "finding reciprocal hits with Python script..."
python find_recip_hits.py > recip_til_AT.txt

echo "done"


