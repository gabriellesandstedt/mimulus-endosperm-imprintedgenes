################################################################################
## filter M. tilingii and M. caespitosa vcfs 
## some filtering methods were modified from: https://evodify.com/gatk-in-non-model-organism/
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 06_filter_vcfs.smk
################################################################################
################################################################################
import os

# assign directories
data_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"
scripts_dir = "/scratch/gds44474/MIMULUS/snps_parents/scripts"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome"

# reference genome: Mimulus IM62 v3
ref = "Mimulus_guttatus_var_IM62.mainGenome.fasta"

# assign samples for individual depth filtering
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# select biallelic SNPs in tilingii vcf 
rule select_biallelic_snps_til:
    input:
        ref=f"{ref_dir}/{ref}",
        vcf=f"{data_dir}/til.vcf"
    output:
        biallelic_vcf=f"{data_dir}/til_biallelic_snps.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -O {output.biallelic_vcf}
        """
        
 # select biallelic SNPs in caes vcf 
rule select_biallelic_snps_caes:
    input:
        ref=f"{ref_dir}/{ref}",
        vcf=f"{data_dir}/caes.vcf"
    output:
        biallelic_vcf=f"{data_dir}/caes_biallelic_snps.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-8.3.0-Java-17.0.4
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -O {output.biallelic_vcf}
        """       
  
 # assign rule to convert vcf to a table to determine quality filtering thresholds with R script
 # problems with java compatability -- had to use different GATK version 
 rule variant_table_til:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/til_biallelic_snps.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics_til.R"
    output:
        table=f"{data_dir}/til_variant.table"
    shell:
        """
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        module load R/4.3.0-foss-2020b
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        R CMD BATCH {input.rscript}
        """
        
   # assign rule to convert vcf to a table to determine quality filtering thresholds with R script
  rule variant_table_caes:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/caes_biallelic_snps.vcf",
        rscript=f"{scripts_dir}/filtering_diagnostics_caes.R"
    output:
        table=f"{data_dir}/caes_variant.table"
    shell:
        """
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        module load R/4.3.0-foss-2020b
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        R CMD BATCH {input.rscript}
        """       
        
# assign rule to filter tilingii vcf 
rule filter_variants_til:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/til_biallelic_snps.vcf"
    output:
        filtered_vcf=f"{data_dir}/til_biallelic_snps_filter.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QUAL < 0 || MQ < 40.00 || SOR > 3.000 || QD < 2.00 || FS > 60.000 || MQRankSum < -12.500 || ReadPosRankSum < -10.000 || ReadPosRankSum > 10.000" \
            --filter-name "my_snp_filter" \
            -O {output.filtered_vcf}
        """
        
# assign rule to filter caespitosa vcf 
rule filter_variants_caes:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/caes_biallelic_snps.vcf"
    output:
        filtered_vcf=f"{data_dir}/caes_biallelic_snps_filter.vcf"
    shell:
        """
        module load gatk/4.1
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QUAL < 0 || MQ < 40.00 || SOR > 3.000 || QD < 2.00 || FS > 60.000 || MQRankSum < -12.500 || ReadPosRankSum < -10.000 || ReadPosRankSum > 10.000" \
            --filter-name "my_snp_filter" \
            -O {output.filtered_vcf}
        """        

# extract passed variants for tilingii vcf 
rule extract_passed_variants_til:
    input:
        filtered_vcf=f"{data_dir}/til_biallelic_snps_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/til_biallelic_snps_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """
        
 # extract passed variants for caespitosa vcf 
rule extract_passed_variants_caes:
    input:
        filtered_vcf=f"{data_dir}/caes_biallelic_snps_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/caes_biallelic_snps_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """   
        
# assign rule to remove heterozygous sites in tilingii vcf 
rule extract_passed_variants_caes:
    input:
        filtered_vcf=f"{data_dir}/caes_biallelic_snps_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/caes_biallelic_snps_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """  

# assign rule to remove heterozygous sites in caespitosa
rule extract_passed_variants_caes:
    input:
        filtered_vcf=f"{data_dir}/caes_biallelic_snps_filter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/caes_biallelic_snps_filterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """  

# assign rule to split vcfs to filter by individual max Depth

# assign rule to filter by depth 

# assign rule to bgzip and index vcfs 

# assign rule to merge zipped vcfs 

# assign rule for final filtering steps of vcf

# assign rule to create final bed files 


        
