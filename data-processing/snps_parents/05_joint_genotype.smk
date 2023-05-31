################################################################################
## joint genotype samples using GATK
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_joint_genotype.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import directory
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome"

# reference genome: Mimulus IM62 v3
ref = "Mimulus_guttatus_var_IM62.mainGenome.fasta"
# define interval list with list of chromosomes
interval_list = "interval.list"
# assign samples
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# define all output files for rule all 
rule all:
    input:
        expand(f"{ref_dir}/{ref}.dict"),
        expand(f"{data_dir}/{{sample}}.g.vcf", sample=samples,
        expand(f"{data_dir}/tilingii.vcf"),
        expand(f"{data_dir}/caes.vcf")

# define rule to index reference with GATK 
rule index_reference:
    input:
        ref = f"{ref_dir}/{ref}"
    output:
        index = f"{ref_dir}/{ref}.dict"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk CreateSequenceDictionary \
            -R {input.ref} \
            -O {output.index}
        """

# define rule to call potential variants for each sample
rule hap_caller:
    input:
        ref = f"{ref_dir}/{ref}",
        bam = f"{data_dir}/{{sample}}_RG_MD_NS_PP_CS.bam",
        intervals = f"{data_dir}/{interval_list}"
    output:
        gvcf = f"{data_dir}/{{sample}}.g.vcf"
    shell:
        """
        module load gatk/4.1
        gatk HaplotypeCaller \
            -I {input.bam} \
            -O {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """

# define rule to combine gvcfs of all samples in a database
rule combine_gvcfs:
    input:
        gvcf=expand(f"{data_dir}/{{sample}}.g.vcf"),
        interval_list=f"{data_dir}/{interval_list}",
    output:
        database=directory("/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/scripts/DB_allsamples")
    shell:
        """
        module load gatk/4.1
        gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.database} \
            --batch-size 158 \
            --reader-threads 6 \
            --sample-name-map {input.map_allsamples} \
            --intervals {input.interval_list}
        """
# define rule to combine gvcfs for retro samples used in the genetic matrix in a database
rule genomicsdb_import_matrix:
    input:
        gvcf=expand(f"{data_dir}/{{sample}}.g.vcf", sample=df2['Sample']),
        interval_list=f"{data_dir}/{interval_list}",
        map_matrix=f"{data_dir}/{sample_map2}"
    output:
        database=directory("/scratch/general/nfs1/u6048240/BOECHERA/GBS_May23/scripts/DB_matrix")
    shell:
        """
        module load gatk/4.1
        gatk GenomicsDBImport \
            --genomicsdb-workspace-path {output.database} \
            --batch-size 41 \
            --reader-threads 6 \
            --sample-name-map {input.map_matrix} \
            --intervals {input.interval_list}
        """

# define rule to joint genotype for all samples at all sites
rule joint_genotype_allsamples:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}"
    output:
        boech_output=f"{data_dir}/boechera_gbs_allsamples.vcf"
    params:
        genomicsdb="gendb://DB_allsamples"

    shell:
        """
        module load gatk/4.1
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {params.genomicsdb} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            --all-sites \
            -O {output.boech_output}
        """
