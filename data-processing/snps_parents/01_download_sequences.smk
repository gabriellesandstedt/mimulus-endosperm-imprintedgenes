################################################################################
## download fastqs from NCBI to call tilingii and caespitosa snps
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 01_download_sequences.smk
################################################################################
################################################################################
import os

# assign data directory
data_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"

# assign samples to be downloaded from NCBI
# "SRR12424410" = SOP12, til
# "SRR3103524" = LVR1, til
# "SRR12424419" = UTC1, caes
# "SRR12424421" = TWN36, caes 
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# define all output files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_1.fastq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_2.fastq.gz", sample=samples)

# define rule to download fastqs from NCBI
# https://www.ncbi.nlm.nih.gov/sra
rule download_fastq:
    output:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    params:
        sra=lambda wildcards: wildcards.sample
    shell:
        """
        echo -e "\\n\\n["$(date)"]\\n Load SRA-Toolkit...\\n\\n"
        module load SRA-Toolkit/2.11.1-centos_linux64
        echo -e "\\n["$(date)"]\\n Download {params.sra} from SRA database \\n"
        prefetch {params.sra} && fastq-dump --split-e --gzip {params.sra} -O {data_dir}
        sleep 60
        """
