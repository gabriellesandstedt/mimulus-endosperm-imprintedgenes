################################################################################
## joint genotype samples using GATK
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_joint_genotype.smk
################################################################################
################################################################################
import os
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/snps_parents_til/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome_til"

# reference genome: Mimulus LVR1 v1 
ref = "Mimulus_tilingii_var_LVR.mainGenome.fasta"
# define interval list with list of chromosomes
interval_list = "intervals.list"
# assign samples
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# define all output files for rule all 
rule all:
    input:
        expand(f"{ref_dir}/Mimulus_tilingii_var_LVR.mainGenome.dict"),
        expand(f"{data_dir}/{{sample}}.g.vcf.gz", sample=samples),
        expand(f"{data_dir}/til.vcf"),
        expand(f"{data_dir}/caes.vcf")

rule index_reference:
    input:
        ref = f"{ref_dir}/{ref}"
    output:
        index = f"{ref_dir}/Mimulus_tilingii_var_LVR.mainGenome.dict"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk CreateSequenceDictionary \
            -R {input.ref} \
            -O {output.index}
        """

# define rule to call potential variants for each sample
rule hap_caller:
    input:
        ref = f"{ref_dir}/{ref}",
        CS_bam = f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam",
        intervals = f"{data_dir}/{interval_list}"
    output:
        gvcf = f"{data_dir}/{{sample}}.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk HaplotypeCaller \
            -I {input.CS_bam} \
            -O {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """

# define rule to combine tilingii gvcfs 
rule combine_gvcfs_til:
    input:
        ref = f"{ref_dir}/{ref}",
        SOP12_gvcf = f"{data_dir}/SRR12424410.g.vcf.gz",
        LVR1_gvcf = f"{data_dir}/SRR3103524.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        til_gvcf = f"{data_dir}/til.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk CombineGVCFs \
            -R {input.ref} \
            --variant {input.SOP12_gvcf} \
            --variant {input.LVR1_gvcf} \
            -L {input.interval_list} \
            -O {output.til_gvcf}
        """
               
# define rule to combine caespitosa gvcfs 
rule combine_gvcfs_caes:
    input:
        ref = f"{ref_dir}/{ref}",
        UTC1_gvcf = f"{data_dir}/SRR12424419.g.vcf.gz",
        TWN36_gvcf = f"{data_dir}/SRR12424421.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        caes_gvcf = f"{data_dir}/caes.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk CombineGVCFs \
            -R {input.ref} \
            --variant {input.UTC1_gvcf} \
            --variant {input.TWN36_gvcf} \
            -L {input.interval_list} \
            -O {output.caes_gvcf}
        """         
               
# define rule to joint genotype for tilingii samples 
rule joint_genotype_til:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}",
        til_gvcf=f"{data_dir}/til.g.vcf.gz"
    output:
        til_vcf=f"{data_dir}/til.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.til_gvcf} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output.til_vcf}
        """
           
# define rule to joint genotype for caespitosa samples 
rule joint_genotype_caes:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}",
        caes_gvcf=f"{data_dir}/caes.g.vcf.gz"
    output:
        caes_vcf=f"{data_dir}/caes.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.caes_gvcf} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output.caes_vcf}
        """
