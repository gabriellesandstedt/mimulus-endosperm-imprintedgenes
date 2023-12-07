# this script is unfinished, just writing my general workflow down # 
wget ftp://ftp.ncbi.nlm.nih.gov//genomes/refseq/plant/Erythranthe_guttata/all_assembly_versions/GCF_000504015.1_Mimgu1_0/GCF_000504015.1_Mimgu1_0_protein.faa.gz
gunzip /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa.gz
makeblastdb -in /scratch/gds44474/MIMULUS/snps_parents_til/data/blastn_db/GCF_000504015.1_Mimgu1_0_protein.faa -dbtype prot
blastx -query shared_34_MEGs.fasta -db /scratch/gds44474/MIMULUS/snps_parents_til/data/blastx_db/GCF_000504015.1_Mimgu1_0_protein.faa  -evalue 1e-10 -max_target_seqs 10 -out shared_34_MEGs_blastx.txt -num_threads 4 -outfmt 6 -max_hsps 1 
