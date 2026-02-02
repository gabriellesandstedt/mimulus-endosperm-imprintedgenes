#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=100:00:00
#SBATCH --nodes=2 #Note: HPCC nodes have 20+ CPU cores each
#SBATCH --ntasks=4 #parallel tasks/processes run on separate CPU cores or nodes (coarse-grained internode/intercore parallelization via MPI; requires OpenMPI/impi) Note: sbatch does not launch the tasks, only requests resources and submits the batch script. ntasks tells Slurm that N parallel tasks will be launched and to allocate resources accordingly. Parallel tasks are launched by the script.
#SBATCH --mem=72G
#SBATCH --job-name callimprinting
#SBATCH --output callimprinting.out
#SBATCH --error callimprinting.error


# SOP is sp A
# LVR is sp B


bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 SOP12xLVR1_1_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -2 LVR1xSOP12_1_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -S diagnostic_snps_tilsp_SOP_LVR_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_lvrsop1

bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 SOP12xLVR1_2_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -2 LVR1xSOP12_2_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -S diagnostic_snps_tilsp_SOP_LVR_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_lvrsop2

bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 SOP12xLVR1_3_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -2 LVR1xSOP12_3_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -S diagnostic_snps_tilsp_SOP_LVR_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_lvrsop3


# UTC is spA
# TWN is spB

bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 UTC1xTWN36_1_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -2 TWN36xUTC1_1_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -S diagnostic_snps_caessp_UTC_TWN_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_twnutc1

bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 UTC1xTWN36_2_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -2 TWN36xUTC1_2_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -S diagnostic_snps_caessp_UTC_TWN_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_twnutc2

bash call_imprinted.sh \
   -s /scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq \
   -1 UTC1xTWN36_3_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \
   -2 TWN36xUTC1_3_STAR_LVR_v1_MD_Split_Q60.noSH.pp.bam \ 
   -S diagnostic_snps_caessp_UTC_TWN_tilref.bed \
   -G Mtilingiivar_LVR_860_v1.1.gene_exons.gtf \
   -R 2 \
   -M 85 \
   -P 50 \
   -o outdir_twnutc3
