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
  
 # assign rule to convert vcf to a table to assess quality filtering thresholds
 # problems with java compatability -- had to use different GATK version 
 rule variant_table_til:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/til_biallelic_snps.vcf",
    output:
        table=f"{data_dir}/til_variant.table"
    shell:
        """
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        """
        
# assign rule to convert vcf to a table to assess quality filtering thresholds
  rule variant_table_caes:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/caes_biallelic_snps.vcf",
    output:
        table=f"{data_dir}/caes_variant.table"
    shell:
        """
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
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
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
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
        module load GATK/4.3.0.0-GCCcore-8.3.0-Java-1.8
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            --filter-expression "QD < 2.0" --filter-name "QD2" \
            --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0" --filter-name "SOR3" \
            --filter-expression "FS > 60.0" --filter-name "FS60" \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
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
rule rem_hets_til:
    input:
        filtered_passed_vcf=f"{data_dir}/til_biallelic_snps_filterPASSED.vcf"
    output:
        nohet_vcf=f"{data_dir}/til_biallelic_snps_filtered_nohets.vcf"
    shell:
        """
        module load BCFtools/1.15.1-GCC-10.2.0
        bcftools view -g ^het {input.filtered_passed_vcf} > {output.nohet_vcf}
        """  
        
# assign rule to remove heterozygous sites in caespitosa vcf 
rule rem_hets_caes:
    input:
        filtered_passed_vcf=f"{data_dir}/caes_biallelic_snps_filterPASSED.vcf"
    output:
        nohet_vcf=f"{data_dir}/caes_biallelic_snps_filtered_nohets.vcf"
    shell:
        """
        module load BCFtools/1.15.1-GCC-10.2.0
        bcftools view -g ^het {input.filtered_passed_vcf} > {output.nohet_vcf}
        """  
        

# assign rule to split til vcf to filter by individual depth
rule split_til_vcf:
    input:
        nohet_vcf=f"{data_dir}/til_biallelic_snps_filtered_nohets.vcf"
    output:
        SOP12_vcf=f"{data_dir}/SOP12_snps.vcf",
        LVR1_vcf=f"{data_dir}/LVR1_snps.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools --remove-indv SRR12424410 --vcf {input.nohet_vcf} --recode --recode-INFO-all --out {output.LVR1_vcf}
        vcftools --remove-indv SRR3103524 --vcf {input.nohet_vcf} --recode --recode-INFO-all --out {output.SOP12_vcf}
        """
        

# assign rule to split caes vcf to filter by individual Depth
rule split_caes_vcf:
    input:
        nohet_vcf=f"{data_dir}/caes_biallelic_snps_filtered_nohets.vcf"
    output:
        UTC1_vcf=f"{data_dir}/UTC1_snps.vcf",
        TWN36_vcf=f"{data_dir}/TWN36_snps.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools --remove-indv SRR12424421 --vcf {input.nohet_vcf} --recode --recode-INFO-all --out {output.UTC1_vcf}
        vcftools --remove-indv SRR12424419 --vcf {input.nohet_vcf} --recode --recode-INFO-all --out {output.TWN36_vcf}
        """ 

# assign rule to filter by depth of tilingii samples 
# mean covg LVR1:  19.8721X ; SD: 307.4309x
# mean covg SOP12: 4.9261X ; SD: 53.9951X
# max DP = 2 x SD + mean 
rule filt_dp_til:
    input:
        SOP12_vcf=f"{data_dir}/SOP12_snps.vcf.recode.vcf",
        LVR1_vcf=f"{data_dir}/LVR1_snps.vcf.recode.vcf"
    output:
        SOP12_dp_vcf=f"{data_dir}/SOP12_snps_dp.vcf",
        LVR1_dp_vcf=f"{data_dir}/LVR1_snps_dp.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools --vcf {input.SOP12_vcf} --maxDP 113 --recode --recode-INFO-all --out {output.SOP12_dp_vcf}
        vcftools --vcf {input.LVR1_vcf} --maxDP 635 --recode --recode-INFO-all --out {output.LVR1_dp_vcf}
        """ 
        
# assign rule to filter by depth of tilingii samples 
# mean covg UTC1:  4.4433X ; SD: 63.0022X
# mean covg TWN36: 5.4595X ; SD: 78.7775X
# max DP = 2 x SD + mean 
rule filt_dp_caes:
    input:
        UTC1_vcf=f"{data_dir}/UTC1_snps.vcf.recode.vcf",
        TWN36_vcf=f"{data_dir}/TWN36_snps.vcf.recode.vcf"
    output:
        UTC1_dp_vcf=f"{data_dir}/UTC1_snps_dp.vcf",
        TWN36_dp_vcf=f"{data_dir}/TWN36_snps_dp.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools --vcf {input.UTC1_vcf} --maxDP 131  --recode --recode-INFO-all --out {output.UTC1_dp_vcf}
        vcftools --vcf {input.TWN36_vcf} --maxDP 163  --recode --recode-INFO-all --out {output.TWN36_dp_vcf}
        """     

# assign rule to bgzip and index vcfs 
rule vcf_to_gzvcf:
    input:
        SOP12_dp_vcf=f"{data_dir}/SOP12_snps_dp.vcf.recode.vcf",
        LVR1_dp_vcf=f"{data_dir}/LVR1_snps_dp.vcf.recode.vcf",
        UTC1_dp_vcf=f"{data_dir}/UTC1_snps_dp.vcf.recode.vcf",
        TWN36_dp_vcf=f"{data_dir}/TWN36_snps_dp.vcf.recode.vcf"
    output:
        SOP12_dp_gzvcf=f"{data_dir}/SOP12_snps_dp.vcf.recode.vcf.gz",
        LVR1_dp_gzvcf=f"{data_dir}/LVR1_snps_dp.vcf.recode.vcf.gz",
        UTC1_dp_gzvcf=f"{data_dir}/UTC1_snps_dp.vcf.recode.vcf.gz",
        TWN36_dp_gzvcf=f"{data_dir}/TWN36_snps_dp.vcf.recode.vcf.gz",
        SOP12_dp_tabix=f"{data_dir}/SOP12_snps_dp.vcf.recode.vcf.gz.tbi",
        LVR1_dp_tabix=f"{data_dir}/LVR1_snps_dp.vcf.recode.vcf.gz.tbi",
        UTC1_dp_tabix=f"{data_dir}/UTC1_snps_dp.vcf.recode.vcf.gz.tbi",
        TWN36_dp_tabix=f"{data_dir}/TWN36_snps_dp.vcf.recode.vcf.gz.tbi"
    shell:
        """
        module load HTSlib/1.15.1-GCC-10.2.0
        bgzip {input.SOP12_dp_vcf}
        bgzip {input.LVR1_dp_vcf}
        bgzip {input.UTC1_dp_vcf}
        bgzip {input.TWN36_dp_vcf}
        tabix -p vcf {output.SOP12_dp_gzvcf}
        tabix -p vcf {output.LVR1_dp_gzvcf}
        tabix -p vcf {output.UTC1_dp_gzvcf}
        tabix -p vcf {output.TWN36_dp_gzvcf}
        """

# assign rule to merge zipped vcfs 
rule merge_vcfs:
    input:
        SOP12_dp_gzvcf=f"{data_dir}/SOP12_snps_dp.vcf.recode.vcf.gz",
        LVR1_dp_gzvcf=f"{data_dir}/LVR1_snps_dp.vcf.recode.vcf.gz",
        UTC1_dp_gzvcf=f"{data_dir}/UTC1_snps_dp.vcf.recode.vcf.gz",
        TWN36_dp_gzvcf=f"{data_dir}/TWN36_snps_dp.vcf.recode.vcf.gz"
    output:
        til_merged_vcf=f"{data_dir}/til_snps_filtered_maxdp.vcf",
        caes_merged_vcf=f"{data_dir}/caes_snps_filtered_maxdp.vcf"
    shell:
        """
        module load  BCFtools/1.15.1-GCC-10.2.0
        bcftools merge {input.SOP12_dp_gzvcf} {input.LVR1_dp_gzvcf}  > {output.til_merged_vcf}
        bcftools merge {input.UTC1_dp_gzvcf} {input.TWN36_dp_gzvcf}  > {output.caes_merged_vcf}
        """   

# assign rule for final filtering steps of tilingii vcf
# min DP 4 
# doesn't allow any missing genotypes with max-missing 0
# minor allele count (mac) of 1
# 1085414 sites remained
rule final_til_vcf:
    input:
        til_merged_vcf=f"{data_dir}/til_snps_filtered_maxdp.vcf"
    output:
        final_til_vcf=f"{data_dir}/til_snps_filtered_DP_mac.vcf",
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools \
         --vcf {input.til_merged_vcf} \
         --remove-indels \
         --minDP 4 \
         --max-missing 0 \
         --mac 1 \
         --recode \
         --recode-INFO-all \
         --out {output.final_til_vcf}
        """  

# assign rule for final filtering steps of tilingii vcf
# min DP 4 
# doesn't allow any missing genotypes with max-missing 0
# minor allele count (mac) of 1
# 322508 sites remained
rule final_caes_vcf:
    input:
        caes_merged_vcf=f"{data_dir}/caes_snps_filtered_maxdp.vcf"
    output:
        final_caes_vcf=f"{data_dir}/caes_snps_filtered_DP_mac.vcf",
    shell:
        """
        module load VCFtools/0.1.16-GCC-10.2.0
        vcftools \
         --vcf {input.caes_merged_vcf} \
         --remove-indels \
         --minDP 4 \
         --max-missing 0 \
         --mac 1 \
         --recode \
         --recode-INFO-all \
         --out {output.final_caes_vcf}
        """ 

# assign rule to create final bed files 


        
