cd import os
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq"
ref_dir = "/scratch/gds44474/MIMULUS/rna_seq_26/til_rnaseq"

# reference genome: Mimulus LVR1 v1 
ref = "Mtilingiivar_LVR_860_v1.0.fa"
# define interval list with list of chromosomes
interval_list = "intervals.list"
# assign samples
samples = ["SRR12424410", "SRR3103524"]

rule all:
    input:
        expand(f"{data_dir}/til.g.vcf.gz")

rule joint_genotype_all:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}",
        all_gvcf=f"{data_dir}/til.g.vcf.gz"
    output:
        all_vcf=f"{data_dir}/til_allsites.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.all_gvcf} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            --all-sites \
            --output {output.all_vcf}
        """
