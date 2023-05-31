################################################################################
## joint genotype samples using GATK
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 05_joint_genotype.smk
################################################################################
################################################################################
import os
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
        gvcf = f"{data_dir}/{{sample}}.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk HaplotypeCaller \
            -I {input.bam} \
            -O {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """

# define rule to combine tilingii gvcfs 
rule combine_gvcfs_til:
    input:
        SOP12_gvcf = f"{data_dir}/SRR12424410.g.vcf.gz",
        LVR1_gvcf = f"{data_dir}/SRR3103524.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        til_gvcf = f"{data_dir}/til.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk CombineGVCFs \
            -R reference.fasta \
            --variant {input.SOP12_gvcf} \
            --variant {input.LVR1_gvcf} \
            -L {input.interval_list} \
            -O {output.til_gvcf}
        """
               
# define rule to combine caespitosa gvcfs 
rule combine_gvcfs_caes:
    input:
        UTC1_gvcf = f"{data_dir}/SRR12424419.g.vcf.gz",
        TWN36_gvcf = f"{data_dir}/SRR12424421.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        caes_gvcf = f"{data_dir}/caes.g.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk CombineGVCFs \
            -R reference.fasta \
            --variant {input.UTC1_gvcf} \
            --variant {input.TWN36_gvcf} \
            -L {input.interval_list} \
            -O {output.caes_gvcf}
        """         
               
# define rule to joint genotype for tilingii samples 
rule joint_genotype_til:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}"
        til_gvcf=f"{data_dir}/til.g.vcf.gz"
    output:
        til_vcf=f"{data_dir}/til.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
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
        intervals=f"{data_dir}/{interval_list}"
        caes_gvcf=f"{data_dir}/caes.g.vcf.gz"
    output:
        caes_vcf=f"{data_dir}/caes.vcf.gz"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk GenotypeGVCFs \
            -R {input.ref} \
            -V {input.caes_gvcf} \
            -L {input.intervals} \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output.caes_vcf}
        """
