# this script is unfinished, just writing my general workflow down # 
wget ftp://ftp.ncbi.nlm.nih.gov//genomes/refseq/plant/Erythranthe_guttata/all_assembly_versions/GCF_000504015.1_Mimgu1_0/GCF_000504015.1_Mimgu1_0_protein.faa.gz
gunzip /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/blastn_db/GCF_000504015.1_Mimgu1_0_protein.faa -dbtype prot
blastx -query shared_34_MEGs.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa  -evalue 1e-10 -max_target_seqs 10 -out shared_34_MEGs_blastx.txt -num_threads 4 -outfmt 6 -max_hsps 1 

wget 
gunzip Araport11_pep_20220914_representative_gene_model.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/Araport11_pep_20220914_representative_gene_model -dbtype prot
blastx -query shared_34_MEGs.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/tair_protein/Araport11_pep_20220914_representative_gene_model  -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out shared_34_MEGs_tair.txt



#https://www.genoscope.cns.fr/brassicanapus/data/
wget https://www.genoscope.cns.fr/brassicanapus/data/Brassica_napus_v4.1.chromosomes.fa.gz
gunzip Brassica_napus.annotation_v5.pep.fa.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus/Brassica_napus.annotation_v5.pep.fa -dbtype prot -out brassicanapus
blastx -query til_197megs_sequences.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus  -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_197megs_sequences_brassica_napus.txt
blastx -query til_197megs_sequences.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/brassicanapus -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -max_hsps 1 -outfmt 6 -out til_197megs_sequences_brassica_napus.txt



# Download and prepare sequences
wget https://www.genoscope.cns.fr/brassicanapus/data/Brassica_napus.annotation_v5.pep.fa.gz
gunzip Brassica_napus.annotation_v5.pep.fa.gz

# Make BLAST databases
makeblastdb -in Brassica_napus.annotation_v5.cds.fa -dbtype nucl -out brassicanapus
makeblastdb -in genes_of_interest.fasta -dbtype nucl -out tilingii

# Perform BLAST searches
blastn -query genes_of_interest.fasta -db brassicanapus -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_brassica.txt
blastn -query Brassica_napus.annotation_v5.cds.fa -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out brassica_to_til.txt

# Run Python script to find reciprocal hits
python find_recip_hits.py


blastx -query genes_of_interest.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/Araport11_pep_20220914_representative_gene_model  -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_AT.txt
tblastn -query /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/Araport11_pep_20220914_representative_gene_model -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out AT_to_til.txt
python find_recip_hits.py > recip_til_AT.txt


makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/BLAST_dbs/Crubella_183_v1.0.cds.fa -dbtype nucl -out capsella_rubella
blastn -query genes_of_interest.fasta -db capsella_rubella -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_caprubella.txt
blastn -query Crubella_183_v1.0.cds.fa -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out caprubella_to_til.txt
python find_recip_hits.py > recip_til_caprubella.txt


makeblastdb -in Alyrata_384_v2.1.cds.fa -dbtype nucl -out arabidopsis_lyrata
blastn -query genes_of_interest.fasta -db arabidopsis_lyrata -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_alyrata.txt
blastn -query Alyrata_384_v2.1.cds.fa -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out alyrata_to_til.txt
python find_recip_hits.py > recip_til_alyrata.txt

makeblastdb -in ZmB73_5a.59_working.cds.fasta -dbtype nucl -out zea_mays
blastn -query genes_of_interest.fasta -db zea_mays -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_zm.txt
blastn -query ZmB73_5a.59_working.cds.fasta -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out zm_to_til.txt
python find_recip_hits.py > zm_til_recip_hits.txt


makeblastdb -in Slycopersicum_225_cds.fa -dbtype nucl -out sol
blastn -query genes_of_interest.fasta -db sol -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out til_to_sol.txt
blastn -query Slycopersicum_225_cds.fa -db tilingii -evalue 1e-10 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out sol_to_til.txt
python find_recip_hits.py > sol_til_recip_hits.txt



