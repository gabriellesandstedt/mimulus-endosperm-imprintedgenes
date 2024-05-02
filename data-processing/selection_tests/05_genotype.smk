cd import os
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/selection_tests/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome_til"

# reference genome: Mimulus LVR1 v1 
ref = "Mimulus_tilingii_var_LVR.mainGenome.fasta"
# define interval list with list of chromosomes
interval_list = "intervals.list"
# assign samples
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421","SRR12424411", "SRR12424412","SRR12424417","SRR12424423", "SRR12424422", "SRR12424418","SRR12424416", "SRR12424413"]


rule all:
    input:
        expand(f"{data_dir}/{{sample}}.g.vcf.gz", sample=samples),
        f"{data_dir}/til.g.vcf.gz",
        f"{data_dir}/caes.g.vcf.gz"

# changed -O to -o for gatk 3.8
# 3.8 also requires -T for calling HaplotypeCaller
# java -jar $EBROOTGATK/GenomeAnalysisTK.jar for 3.8
# gatk 4 only requires "gatk HaplotypeCaller"

rule hap_caller:
    input:
        ref = f"{ref_dir}/{ref}",
        CS_bam = f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam",
        intervals = f"{data_dir}/{interval_list}"
    output:
        gvcf = f"{data_dir}/{{sample}}.g.vcf.gz"
    shell:
        """
        module load GATK/3.8-1-Java-1.8.0_144
        gatk HaplotypeCaller \
            -I {input.CS_bam} \
            -o {output.gvcf} \
            -R {input.ref} \
            -L {input.intervals} \
            -ERC GVCF
        """

# define rule to combine tilingii gvcfs 
# I renamed the vcfs from their SRA IDs to their population & line names 
rule combine_gvcfs_til:
    input:
        ref = f"{ref_dir}/{ref}",
        SOP12_gvcf = f"{data_dir}/SOP12.g.vcf.gz",
        LVR1_gvcf = f"{data_dir}/LVR.g.vcf.gz",
        ICE10_gvcf = f"{data_dir}/ICE10.g.vcf.gz",
        SAB1_gvcf = f"{data_dir}/SAB1.g.vcf.gz",
        SAB19_gvcf = f"{data_dir}/SAB19.g.vcf.gz",
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
            --variant {input.ICE10_gvcf} \
            --variant {input.SAB1_gvcf} \
            --variant {input.SAB19_gvcf} \
            -L {input.interval_list} \
            -O {output.til_gvcf}
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
            --all-sites \
            -O {output.til_vcf}
        """


# define rule to combine caes gvcfs 
# I renamed the vcfs from their SRA IDs to their population & line names 
rule combine_gvcfs_caes:
    input:
        ref = f"{ref_dir}/{ref}",
        UTC1_gvcf = f"{data_dir}/UTC1.g.vcf.gz",
        UTC2_gvcf = f"{data_dir}/UTC2.g.vcf.gz",
        GAB1_gvcf = f"{data_dir}/GAB1.g.vcf.gz",
        GAB2_gvcf = f"{data_dir}/GAB2.g.vcf.gz",
        PAG2_gvcf = f"{data_dir}/PAG2.g.vcf.gz",
        KCK1_gvcf = f"{data_dir}/KCK1.g.vcf.gz",
        TWN36_gvcf = f"{data_dir}/TWN36.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        caes_gvcf = f"{data_dir}/caes.g.vcf.gz"  
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk CombineGVCFs \
            -R {input.ref} \
            --variant {input.UTC1_gvcf} \
            --variant {input.UTC2_gvcf} \
            --variant {input.GAB1_gvcf} \
            --variant {input.GAB2_gvcf} \
            --variant {input.PAG2_gvcf} \
            --variant {input.KCK1_gvcf} \
            --variant {input.TWN36_gvcf} \
            -L {input.interval_list} \
            -O {output.caes_gvcf}
        """


# define rule to joint genotype for caes samples 
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
            --all-sites \
            -O {output.caes_vcf}
        """

# define rule to combine all gvcfs 
# I renamed the vcfs from their SRA IDs to their population & line names 
rule combine_gvcfs_caes:
    input:
        ref = f"{ref_dir}/{ref}",
        UTC1_gvcf = f"{data_dir}/UTC1.g.vcf.gz",
        UTC2_gvcf = f"{data_dir}/UTC2.g.vcf.gz",
        GAB1_gvcf = f"{data_dir}/GAB1.g.vcf.gz",
        GAB2_gvcf = f"{data_dir}/GAB2.g.vcf.gz",
        PAG2_gvcf = f"{data_dir}/PAG2.g.vcf.gz",
        KCK1_gvcf = f"{data_dir}/KCK1.g.vcf.gz",
        TWN36_gvcf = f"{data_dir}/TWN36.g.vcf.gz",
        SOP12_gvcf = f"{data_dir}/SOP12.g.vcf.gz",
        LVR1_gvcf = f"{data_dir}/LVR.g.vcf.gz",
        ICE10_gvcf = f"{data_dir}/ICE10.g.vcf.gz",
        SAB1_gvcf = f"{data_dir}/SAB1.g.vcf.gz",
        SAB19_gvcf = f"{data_dir}/SAB19.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        all_gvcf = f"{data_dir}/til_caes.g.vcf.gz"  
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk CombineGVCFs \
            -R {input.ref} \
            --variant {input.UTC1_gvcf} \
            --variant {input.UTC2_gvcf} \
            --variant {input.GAB1_gvcf} \
            --variant {input.GAB2_gvcf} \
            --variant {input.PAG2_gvcf} \
            --variant {input.KCK1_gvcf} \
            --variant {input.TWN36_gvcf} \
            --variant {input.SOP12_gvcf} \
            --variant {input.LVR1_gvcf} \
            --variant {input.ICE10_gvcf} \
            --variant {input.SAB1_gvcf} \
            --variant {input.SAB19_gvcf} \
            -L {input.interval_list} \
            -O {output.all_gvcf}
        """



rule combine_gvcfs_caes:
    input:
        ref = f"{ref_dir}/{ref}",
        UTC1_gvcf = f"{data_dir}/UTC1_gatk3.g.vcf.gz",
        UTC2_gvcf = f"{data_dir}/UTC2_gatk3.g.vcf.gz",
        GAB1_gvcf = f"{data_dir}/GAB1_gatk3.g.vcf.gz",
        GAB2_gvcf = f"{data_dir}/GAB2_gatk3.g.vcf.gz",
        PAG2_gvcf = f"{data_dir}/PAG2_gatk3.g.vcf.gz",
        KCK1_gvcf = f"{data_dir}/KCK1_gatk3.g.vcf.gz",
        TWN36_gvcf = f"{data_dir}/TWN36_gatk3.g.vcf.gz",
        SOP12_gvcf = f"{data_dir}/SOP12_gatk3.g.vcf.gz",
        LVR1_gvcf = f"{data_dir}/LVR1_gatk3.g.vcf.gz",
        ICE10_gvcf = f"{data_dir}/ICE10_gatk3.g.vcf.gz",
        SAB1_gvcf = f"{data_dir}/SAB1_gatk3.g.vcf.gz",
        SAB19_gvcf = f"{data_dir}/SAB19_gatk3.g.vcf.gz",
        interval_list = f"{data_dir}/{interval_list}"
    output:
        all_gvcf = f"{data_dir}/til_caes_gatk3.g.vcf.gz"
    shell:
        """
        module load GATK/3.8-1-Java-1.8.0_144
        java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
            -T CombineGVCFs \
            -R {input.ref} \
            --variant {input.UTC1_gvcf} \
            --variant {input.UTC2_gvcf} \
            --variant {input.GAB1_gvcf} \
            --variant {input.GAB2_gvcf} \
            --variant {input.PAG2_gvcf} \
            --variant {input.KCK1_gvcf} \
            --variant {input.TWN36_gvcf} \
            --variant {input.SOP12_gvcf} \
            --variant {input.LVR1_gvcf} \
            --variant {input.ICE10_gvcf} \
            --variant {input.SAB1_gvcf} \
            --variant {input.SAB19_gvcf} \
            -L {input.interval_list} \
            -o {output.all_gvcf}
        """


rule joint_genotype_all:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}",
        all_gvcf=f"{data_dir}/til_caes.g.vcf.gz"
    output:
        all_vcf=f"{data_dir}/til_caes.vcf"
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


rule joint_genotype_all:
    input:
        ref=f"{ref_dir}/{ref}",
        intervals=f"{data_dir}/{interval_list}",
        all_gvcf=f"{data_dir}/til_caes_gatk3.g.vcf.gz"
    output:
        all_vcf=f"{data_dir}/til_caes_gatk3.vcf"
    shell:
        """
        module load GATK/3.8-1-Java-1.8.0_144
            java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
            -T GenotypeGVCFs \
            -R {input.ref} \
            -V {input.all_gvcf} \
            -L {input.intervals} \
            --standard_min_confidence_threshold_for_calling 30 \
            --includeNonVariantSites \
            -o {output.all_vcf}
        """

