################################################################################
## Filter bams with picard, GATK, and samtools 
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 03_filter_bams.smk
################################################################################
################################################################################
import os 

# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/rna_seq_26/caes_rnaseq/star_pass"
ref_dir = "/scratch/gds44474/MIMULUS/rna_seq_26/caes_rnaseq"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]

# assign reference genome
ref = "Mcaespitosavar_TWN36_992_v1.1.fa"

# define all output files to rule all
rule all:
    input:
        f"{ref_dir}/{ref}.dict",
        expand(f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split_Q60.bam", sample=samples)

# create index file for reference genome
rule index_reference:
    input:
        fa2=f"{ref_dir}/{ref}"
    output:
        ref_index=f"{ref_dir}/{ref}.dict"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-12.3.0-Java-17
        gatk CreateSequenceDictionary \
            -R {input.fa2} \
            -O {output.ref_index}
        """
        
# define rule to remove unmapped (4), QC-fail (512), secondary (256), supplementary (2048) and sort bam files 
# samtools v 1.16: https://github.com/samtools/samtools
rule sort_and_index_bam:
    input:
        bam=f"{star_pass2_dir}/{{sample}}.Aligned.sortedByCoord.out.bam"
    output:
        sorted_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1.bam",
        sorted_bai=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1.bam.bai"
    shell:
        """
        module load SAMtools/1.21-GCC-13.3.0
        samtools view -hu -F 2820 {input.bam} | samtools sort -O bam -o {output.sorted_bam} -T {output.sorted_bam}.tmp -
        samtools index {output.sorted_bam}
        """


# define rule to mark duplicates
# picard v 2.27: https://broadinstitute.github.io/picard/
rule mark_duplicates:
    input:
        sorted_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1.bam"
    output:
        MD_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD.bam",
        metrics_file=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD.txt"
    shell:
        """
        module load picard/3.3.0-Java-17
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT={input.sorted_bam} OUTPUT={output.MD_bam} METRICS_FILE={output.metrics_file} REMOVE_DUPLICATES=FALSE CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT ASSUME_SORT_ORDER=coordinate
        """

# define rule to trim intronic spanning reads and automatically transform default star mapq value for uniquely mapping reads of 225 to 60
# GATK v 4.4: https://gatk.broadinstitute.org/hc/en-us
rule split_trim_reads:
    input:
        MD_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD.bam",
        fa2 = f"{ref_dir}/{ref}"
    output:
        split_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split.bam"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-12.3.0-Java-17
        gatk SplitNCigarReads -I {input.MD_bam} -O {output.split_bam} -R {input.fa2}
        """

# define rule to filter reads with a MAPQ score below 60
# samtools v 1.16: https://github.com/samtools/samtools
rule filter_unique_reads:
    input:
        split_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split.bam"
    output:
        filtered_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split_Q60.bam"
    shell:
        """
        module load SAMtools/1.21-GCC-13.3.0
        samtools view -hu -q 60 {input.split_bam} | samtools sort -O bam -o {output.filtered_bam} -T {output.filtered_bam}.tmp -
        samtools index {output.filtered_bam}
        """

# define rule to check details of bam files 
# samtools v 1.16: https://github.com/samtools/samtools
rule flagstat_quality_check:
    input:
        filtered_bam=f"{star_pass2_dir}/{{sample}}_STAR_TWN_v1_MD_Split_Q60.bam"
    shell:
        """
        module load SAMtools/1.21-GCC-13.3.0
        samtools flagstat {input.filtered_bam}
        """
