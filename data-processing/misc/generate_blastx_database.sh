# this script is unfinished, just writing my general workflow down # 
wget ftp://ftp.ncbi.nlm.nih.gov//genomes/refseq/plant/Erythranthe_guttata/all_assembly_versions/GCF_000504015.1_Mimgu1_0/GCF_000504015.1_Mimgu1_0_protein.faa.gz
gunzip /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/blastn_db/GCF_000504015.1_Mimgu1_0_protein.faa -dbtype prot
blastx -query shared_34_MEGs.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa  -evalue 1e-10 -max_target_seqs 10 -out shared_34_MEGs_blastx.txt -num_threads 4 -outfmt 6 -max_hsps 1 


gunzip Araport11_pep_20220914_representative_gene_model.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model -dbtype prot
blastx -query shared_34_MEGs.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model  -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out shared_34_MEGs_tair.txt



#https://www.genoscope.cns.fr/brassicanapus/data/
wget https://www.genoscope.cns.fr/brassicanapus/data/Brassica_napus.annotation_v5.pep.fa.gz
gunzip Brassica_napus.annotation_v5.pep.fa.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus/Brassica_napus.annotation_v5.pep.fa -dbtype prot -out brassicanapus
blastx -query til_197megs_sequences.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus  -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_197megs_sequences_brassica_napus.txt
blastx -query til_197megs_sequences.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -max_hsps 1 -outfmt 6 -out til_197megs_sequences_brassica_napus.txt



# Download and prepare sequences
wget https://www.genoscope.cns.fr/brassicanapus/data/Brassica_napus.annotation_v5.pep.fa.gz
gunzip Brassica_napus.annotation_v5.pep.fa.gz

# Make BLAST databases
makeblastdb -in Brassica_napus.annotation_v5.pep.fa -dbtype prot -out brassicanapus
makeblastdb -in genes_of_interest.fasta -dbtype nucl -out other_species

# Perform BLAST searches
blastx -query genes_of_interest.fasta -db brassicanapus -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_brassica.txt
tblastn -query Brassica_napus.annotation_v5.pep.fa -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out brassica_to_til.txt

# Run Python script to find reciprocal hits
python find_reciprocal_hits.py
