################################################################################
## rename samples to species 
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_feature_counts.smk
################################################################################
################################################################################
import os

# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq"
ref_dir="/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq"

# assign gff
gff = "Mtilingiivar_LVR_860_v1.1.gene_exons.gff3"

def get_output_filenames():
    prefix_mapping = {
        "13_S17": "LVR1xSOP12_1",
        "41_S24": "LVR1xSOP12_2",
        "50_S30": "LVR1xSOP12_3",
        "15_S7": "SOP12xLVR1_1",
        "39_S23": "SOP12xLVR1_2",
        "46_S26": "SOP12xLVR1_3",
        "35_S10": "TWN36xUTC1_1",
        "52_S15": "TWN36xUTC1_2",
        "45_S14": "TWN36xUTC1_3",
        "31_S8": "UTC1xTWN36_1",
        "33_S22": "UTC1xTWN36_2",
        "48_S28": "UTC1xTWN36_3"
    }
    for prefix in ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28"]:
        new_prefix = prefix_mapping[prefix]
        yield f"{star_pass2_dir}/{new_prefix}_STAR_LVR_v1_MD_Split_Q60.bam"

rule rename_files:
    input:
        expand("{prefix}_STAR_LVR_v1_MD_Split_Q60.bam", prefix=["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28"], star_pass2_dir=star_pass2_dir)
    output:
        output_files = list(get_output_filenames())
    run:
        for old_file, new_file in zip(input, output.output_files):
            shell("mv {old_file} {new_file}")
