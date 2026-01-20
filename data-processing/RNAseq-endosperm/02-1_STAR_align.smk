################################################################################
## index and align with STAR
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 02-1_STAR_align.smk
################################################################################
################################################################################
import os


### this v2 script does not deal with masking genome ####
# assign directories
STAR_genome_dir="/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq/STARgenome"
data_dir="/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq"
star_pass_dir="/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq/star_pass"

# assign genome files
ref = "Mtilingiivar_LVR_860_v1.0.fa"
gtf = "Mtilingiivar_LVR_860_v1.1.gene_exons.gtf"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]

# define all output files to rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=samples)

rule create_star_index:
    input:
        fa = f"{data_dir}/{ref}",
        gtf = f"{data_dir}/{gtf}"
    output:
        genome_sentinel = f"{STAR_genome_dir}/Genome"
    threads: 4
    shell:
        """
        ml STAR/2.7.10b-GCC-11.3.0
        mkdir -p {STAR_genome_dir}
        STAR --runThreadN {threads} --runMode genomeGenerate \
             --genomeDir {STAR_genome_dir} --genomeFastaFiles {input.fa} \
             --sjdbGTFfile {input.gtf} \
             --genomeSAindexNbases 13 --sjdbOverhang 149
        """

rule star_alignment_2pass_basic:
    input:
        genome_sentinel = f"{STAR_genome_dir}/Genome",
        rd1=f"{data_dir}/{{sample}}_R1_trim.fastq.gz",
        rd2=f"{data_dir}/{{sample}}_R2_trim.fastq.gz"
    output:
        bam=f"{data_dir}/{{sample}}.Aligned.sortedByCoord.out.bam",
        log_final=f"{data_dir}/{{sample}}.Log.final.out",
        sj=f"{data_dir}/{{sample}}.SJ.out.tab"
    params:
        prefix = f"{data_dir}/{{sample}}."
    threads: 4
    shell:
        """
        ml STAR/2.7.10b-GCC-11.3.0
        STAR --runThreadN {threads} \
             --genomeDir {STAR_genome_dir} \
             --twopassMode Basic \
             --readFilesCommand zcat \
             --readFilesIn {input.rd1} {input.rd2} \
             --alignIntronMin 20 \
             --alignIntronMax 10000 \
             --outFilterMismatchNoverReadLmax 0.05 \
             --outSAMmapqUnique 60 \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes All \
             --outSAMattrRGline ID:{wildcards.sample} LB:LVR.v1 SM:{wildcards.sample} PL:Illumina
        """
