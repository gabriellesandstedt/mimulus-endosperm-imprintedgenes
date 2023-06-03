################################################################################
## Index reference genome, align trimmed fastq files, remove reads with a map quality score < 29
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 03_align.smk
################################################################################
################################################################################
import os
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome"

# reference genome: Mimulus IM62 v3
ref = "Mimulus_guttatus_var_IM62.mainGenome.fasta"

# assign samples to be aligned
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# define all output files from script
rule all:
    input:
        expand(f"{ref_dir}/{ref}.{{suffix}}", suffix=["bwt", "amb", "ann", "pac"], strict=False),
        expand(f"{data_dir}/{{sample}}.sam", sample=samples),
        expand(f"{data_dir}/{{sample}}.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_sorted.bam", sample=samples) 

# define rule to index reference genome using bwa index
# BWA v 0.7.17 : https://bio-bwa.sourceforge.net
rule bwa_index:
    input:
        ref = f"{ref_dir}/{ref}"
    output:
        bwt = "{ref}.bwt",
        amb = "{ref}.amb",
        ann = "{ref}.ann",
        pac = "{ref}.pac"
    shell:
        """
        module load BWA/0.7.17-GCC-10.3.0
        bwa index -a bwtsw {input.ref}
        echo -e "\n["$(date)"]\n done ...\n"
        """

# define rule to align fastqs to reference genome
# BWA v 0.7.17: https://bio-bwa.sourceforge.net
rule bwa_mem:
    input:
        ref = f"{ref_dir}/{ref}",
        trim_fq = f"{data_dir}/{{sample}}_trim.fq.gz"
    output:
        sam = f"{data_dir}/{{sample}}.sam"
    shell:
        """
        module load BWA/0.7.17-GCC-10.3.0
        echo -e "\\n["$(date)"]\\n run BWA mem..\\n"
        bwa mem {input.ref} {input.trim_fq} > {output.sam}
        """

# define rule to convert raw alignments to bam files and removing reads < map quality 29
# Samtools v 1.16 : https://github.com/samtools/samtools
rule sam_bam_q29:
    input:
        sam = f"{data_dir}/{{sample}}.sam"
    output:
        bam = f"{data_dir}/{{sample}}.bam"
    shell:
        """
        module load samtools/1.16
        samtools view -q 29 -b {input.sam} > {output.bam}
        echo -e "\\n["$(date)"]\\n reads with map quality <29 are removed..\\n"
        """

# define rule to sort bam file
# Samtools v 1.16 : https://github.com/samtools/samtools
rule sort_bam:
    input:
        bam = f"{data_dir}/{{sample}}.bam"
    output:
        sorted_bam = f"{data_dir}/{{sample}}_sorted.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort {input.bam} -o {output.sorted_bam}
        samtools index {input.sorted_bam}
        echo -e "\\n["$(date)"]\\n bam file is sorted ..\\n"
        """

