################################################################################
## extract caespitosa and tilingii feature counts
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_feature_counts.smk
################################################################################
################################################################################
import os

# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"
scripts_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm/scripts"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome"

# assign files
gff = "MguttatusvarIM62v3.1.primaryTrs.gff3"

samp_nos = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28"]
samp_names = ["LVR1xSOP12_1", "LVR1xSOP12_2", "LVR1xSOP12_3", "SOP12xLVR1_1", "SOP12xLVR1_2", "SOP12xLVR1_3", "TWN36xUTC1_1", "TWN36xUTC1_2", "TWN36xUTC1_3", "UTC1xTWN36_1", "UTC1xTWN36_2", "UTC1xTWN36_3"]

rule all:
    input:
        expand(f"{star_pass2_dir}/{{samp_no}}_STAR_IM62_v3_MD_Split_Q60_NS.bam", samp_no=samp_nos)
        expand(f"{star_pass2_dir}/{{samp_name}}_STAR_IM62_v3_MD_Split_Q60_NS.bam", samp_name=samp_names)
rule assign:
    output:
        samp_no=samp_nos
    shell:
        """
        print {output}
        """
rule rename_files:
    input:
        f"{star_pass2_dir}/{{samp_no}}_STAR_IM62_v3_MD_Split_Q60_NS.bam", samp_no=sampe_nos
    output:
        f"{star_pass2_dir}/{samp_name}_STAR_IM62_v3_MD_Split_Q60_NS.bam", samp_name=samp_names
    shell:
        "cp {input} {output}"

rule feature_counts:
    input:
        rscript=f"{scripts_dir}/feature_counts.r",
        gff=f"{ref_dir}/{gff}"
    output:
        caes_counts=f"{star_pass2_dir}/endosperm_counts_caes_object",
        til_counts=f"{star_pass2_dir}/endosperm_counts_til_object"
    shell:
        """
        module load R/4.2.1-foss-2020b
        R CMD BATCH {input.rscript} \
        --annotation="{input.gff_file}" \
        --output="{star_pass2_dir}/"
        """

