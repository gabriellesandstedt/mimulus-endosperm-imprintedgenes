data_dir = "/scratch/gds44474/MIMULUS/selection_tests/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome_til"

# reference genome: Mimulus LVR1 v1 
ref = "Mimulus_tilingii_var_LVR.mainGenome.fasta"


# select biallelic SNPs
rule select_biallelic_snps:
    input:
        ref=f"{ref_dir}/{ref}",
        vcf=f"{data_dir}/til_caes.vcf"
    output:
        biallelic_vcf=f"{data_dir}/til_caes_biallelic_snps.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            -O {output.biallelic_vcf}
        """

# select invariant sites
rule select_invariant_sites:
    input:
        ref=f"{ref_dir}/{ref}",
        vcf=f"{data_dir}/til_caes.vcf"
    output:
        invariant_vcf=f"{data_dir}/til_caes_invariant.vcf"
    shell:
        """
        module load  GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include NO_VARIATION \
            -O {output.invariant_vcf}
        """


rule variant_table:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/til_caes_biallelic_snps.vcf",
        rscript=f"{data_dir}/filtering_diagnostics.R"
    output:
        table=f"{data_dir}/til_caes_variant.table"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.biallelic_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """

rule invariant_table:
    input:
        ref=f"{ref_dir}/{ref}",
        invariant_vcf=f"{data_dir}/til_caes_invariant_geno_called.vcf",
        rscript=f"{data_dir}/filtering_diagnostics.R"
    output:
        table=f"{data_dir}/til_caes_invariant.table"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        module load R
        gatk VariantsToTable \
            -R {input.ref} \
            -V {input.invariant_vcf} \
            -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
            -O {output.table}
        
        Rscript {input.rscript}
        """


# rule to filter variants:
rule filter_variants:
    input:
        ref=f"{ref_dir}/{ref}",
        biallelic_vcf=f"{data_dir}/til_caes_biallelic_snps.vcf"
    output:
        filtered_vcf=f"{data_dir}/til_caes_biallelic_snps_qualfilter.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
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

# extract passed variants
rule extract_passed_variants:
    input:
        filtered_vcf=f"{data_dir}/til_caes_biallelic_snps_qualfilter.vcf"
    output:
        filtered_passed_vcf=f"{data_dir}/til_caes_biallelic_snps_qualfilterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_vcf} > {output.filtered_passed_vcf}
        """

# rule to filter invariants:
rule filter_invariants:
    input:
        ref=f"{ref_dir}/{ref}",
        invcf=f"{data_dir}/til_caes_invariant_geno_called.vcf"
    output:
        filtered_invcf=f"{data_dir}/til_caes_invariant_qualfilter.vcf"
    shell:
        """
        module load GATK/4.4.0.0-GCCcore-11.3.0-Java-17
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.invcf} \
            --filter-expression "MQ < 40.0" --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            -O {output.filtered_invcf}
        """

# extract passed invariants
rule extract_passed_invariants:
    input:
        filtered_invcf=f"{data_dir}/til_caes_invariant_qualfilter.vcf"
    output:
        filtered_passed_invcf=f"{data_dir}/til_caes_invariant_qualfilterPASSED.vcf"
    shell:
        """
        grep -E '^#|PASS' {input.filtered_invcf} > {output.filtered_passed_invcf}
        """

samples = ["SRR12424410", "SRR12424411", "SRR12424412", "SRR12424413", "SRR12424416", "SRR12424417", "SRR12424418", "SRR12424419", "SRR12424421", "SRR12424422", "SRR12424423", "SRR3103524"]

# split snp vcf for individual depth filtering
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_snps.vcf", sample=samples)

rule split_vcf:
    input:
        snp_vcf=f"{data_dir}/til_caes_biallelic_snps_qualfilterPASSED.vcf"
    output:
        ind_vcf=f"{data_dir}/{{sample}}_snps.vcf",
    params:
        sample=lambda wildcards: wildcards.sample
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools --indv {params.sample} --vcf {input.snp_vcf} --recode --recode-INFO-all --out {output.ind_vcf}
        """


# split invariant vcf for individual depth filtering
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_invar.vcf", sample=samples)

rule split_vcf:
    input:
        invar_vcf=f"{data_dir}/til_caes_invariant_qualfilterPASSED.vcf"
    output:
        ind_vcf=f"{data_dir}/{{sample}}_invar.vcf",
    params:
        sample=lambda wildcards: wildcards.sample
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools --indv {params.sample} --vcf {input.invar_vcf} --recode --recode-INFO-all --out {output.ind_vcf}
        """


