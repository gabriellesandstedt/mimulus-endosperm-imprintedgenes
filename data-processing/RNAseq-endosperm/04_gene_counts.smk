################################################################################
## 
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 04_gene_counts.smk
################################################################################
################################################################################
import os

# assign directories
ref_dir="/scratch/gds44474/MIMULUS/ref_genome"
star_pass2_dir="/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"

# assign genome files
gff = "MguttatusvarIM62v3.1.primaryTrs.gff3"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]


rule all:
    input:
        expand(f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60_NS.bam" , sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_HTSeq_gene_counts.txt", sample=samples),
        f"{star_pass2_dir}/merged_counts.txt",
        f"{star_pass2_dir}/merged_gene_counts.txt" 
        
rule name_sort_bam:
    input:
        filtered_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60.bam"
    output:
        ns_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60_NS.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        samtools sort {input.filtered_bam} -n -O bam -o {output.ns_bam}
        """
        
rule gene_counts:        
    input:
        ns_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60_NS.bam",
        gff_file=f"{ref_dir}/{gff}"
    output:
        gene_counts=f"{star_pass2_dir}/{{sample}}_HTSeq_gene_counts.txt"
    shell:
        """
        module load HTSeq/0.13.5-foss-2019b-Python-3.7.4
        htseq-count --format bam --stranded no --type gene --idattr Name --nonunique none {input.ns_bam} {input.gff_file} > {output.gene_counts}
        """
        
rule merge_gene_counts:
    input:
        gene_counts=expand(f"{star_pass2_dir}/{{sample}}_HTSeq_gene_counts.txt", sample=samples)
    output:
        merged_counts=f"{star_pass2_dir}/merged_counts.txt"
    shell:
        """
        echo "Column2 {input.gene_counts}" > {output.merged_counts}
        paste -d'\t' {input.gene_counts} | cut -f2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38 >> {output.merged_counts}
        """
        
rule merge_columns_again:
    input:
        gene_col=f"{star_pass2_dir}/15_S7_HTSeq_gene_counts.txt",
        merged_counts=f"{star_pass2_dir}/merged_counts.txt"
    output:
        merged_counts_gene=f"{star_pass2_dir}/merged_gene_counts.txt" 
    shell:
        """
        paste -d'\t' <(cut -f1 {input.gene_col}) {input.merged_counts} > {output.merged_counts_gene}
        """
