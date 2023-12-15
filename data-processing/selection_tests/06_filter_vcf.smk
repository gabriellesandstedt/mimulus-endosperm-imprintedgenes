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
        ref=f"{ref_dir}/{ref_genome}",
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




