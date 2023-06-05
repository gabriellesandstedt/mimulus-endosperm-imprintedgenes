################################################################################
## Filter bams with picard, GATK, and samtools 
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 03_filter_bams.smk
################################################################################
################################################################################
import os 

# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"
repeat_masker_dir = "/scratch/gds44474/MIMULUS/ref_genome/RepeatMasker"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]

# assign reference genome
masked_ref2 = "Mimulus_guttatus_var_IM62.mainGenome.masked.fasta"

# define all output files to rule all
rule all:
    input:
        f"{repeat_masker_dir}/{masked_ref2}.dict",
        expand(f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60.bam", sample=samples)

# create index file for reference genome
rule index_reference:
    input:
        masked_fa2=f"{repeat_masker_dir}/{masked_ref2}"
    output:
        ref_index=f"{repeat_masker_dir}/{masked_ref2}.dict"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk CreateSequenceDictionary \
            -R {input.masked_ref2} \
            -O {output.ref_index}
        """
        
# define rule to remove secondary alignments (-F 524) and sort bam files 
# samtools v 1.16: https://github.com/samtools/samtools
rule sort_and_index_bam:
    input:
        bam=f"{star_pass2_dir}/{{sample}}.bamAligned.sortedByCoord.out.bam"
    output:
        out_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3.bam"
        star_pass2_dir=star_pass2_dir

    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        samtools view -hu -F 524 {input.bam} | samtools sort -O {output.star_pass2_dir} -o {output.bam} -T {output.out_bam}.tmp -
        samtools index {output.bam}
        """

# define rule to mark and remove duplicates
# picard v 2.27: https://broadinstitute.github.io/picard/
rule mark_duplicates:
    input:
        bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3.bam"
    output:
        MD_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD.bam",
        metrics_file=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD.txt"
    shell:
        """
        module load picard/2.27.4-Java-13.0.2
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.MD_bam} METRICS_FILE={output.metrics_file} REMOVE_DUPLICATES=TRUE CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT ASSUME_SORT_ORDER=coordinate
        """

# define rule to trim intronic spanning reads and automatically transform default star mapq value for uniquely mapping reads of 225 to 60
# GATK v 4.4: https://gatk.broadinstitute.org/hc/en-us
rule split_trim_reads:
    input:
        MD_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD.bam"
        masked_fa2 = f"{repeatmasker_dir}/{masked_ref2}"
    output:
        split_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split.bam"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk SplitNCigarReads -I {input.MD_bam} -O {output.split_bam} -R {input.masked_fa2}
        """

# define rule to filter reads with a MAPQ score below 60
# samtools v 1.16: https://github.com/samtools/samtools
rule filter_unique_reads:
    input:
        split_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split.bam"
    output:
        filtered_bam=f"{star_pass2_dir}/{{sample}}_STAR_IM62_v3_MD_Split_Q60.bam"
    shell:
        """
        samtools view -hu -q 60 {input.split_bam} | samtools sort -O bam -o {output.filtered_bam} -T {output.filtered_bam}.tmp -
        samtools index {output.filtered_bam}
        """

# define rule to check details of bam files 
# samtools v 1.16: https://github.com/samtools/samtools
rule flagstat_quality_check:
    input:
        filtered_bam=f"{star_pass2_dir}/{{sample}}.STAR.IM62_v3.DupMark.Split.Q60.bam"
    shell:
        """
        samtools flagstat {input.bam}
        """
