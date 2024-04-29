################################################################################
## download tilingii and caespitosa fastqs from NCBI to perform selection tests
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: snakemake/7.22.0-foss-2022a
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 01_download_sequences.smk
################################################################################
################################################################################
import os

# assign data directory
data_dir = "/scratch/gds44474/MIMULUS/selection_tests/data"

# assign samples to be downloaded from NCBI
# tilingii samples:
# "SRR3103524" = LVR1
# "SRR12424410" = SOP12
# "SRR12424411" = SAB19
# "SRR12424412" = SAB1
# "SRR12424417" = ICE10 

# caespitosa samples: 
# "SRR12424419" = UTC1
# "SRR12424421" = TWN36
# "SRR12424423" = GAB1
# "SRR12424422" = GAB2
# "SRR12424418" = UTC2
# "SRR12424416" = KCK1
# "SRR12424413" = PAG2 
# "SRX142379" = AHQT
# "SRX142377" = SLP

#samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421","SRR12424411", "SRR12424412","SRR12424417","SRR12424423", "SRR12424422", "SRR12424418","SRR12424416", "SRR12424413","SRX142379","SRX142377"]

samples = ["SRX142379","SRX142377"]

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
        module load SRA-Toolkit/3.0.3-gompi-2022a
        echo -e "\\n["$(date)"]\\n Download {params.sra} from SRA database \\n"
        prefetch {params.sra} && fastq-dump --split-3 --gzip {params.sra} -O {data_dir}
        sleep 60
        """
