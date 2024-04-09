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
        mv {output.ind_vcf}.recode.vcf {output.ind_vcf}
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
        mv {output.ind_vcf}.recode.vcf {output.ind_vcf}
        """

samples = ["SRR12424410", "SRR12424411", "SRR12424412", "SRR12424413", "SRR12424416", "SRR12424417", "SRR12424418", "SRR12424419", "SRR12424421", "SRR12424422", "SRR12424423", "SRR3103524"]
maxDP = ["113","66","123","133","199","105","115","130","163","149","130","635"]
#maxDP_lvr1 = ["127","68","137","134","213","123", "124", "138", "183", "164", "137","671"]
#maxDP_im62 = ["113","66","123","133","199","105","115","130","163","149","130","635"]
# filter for individual depth of snp files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf", sample=samples)

rule filt_dp_snp_samples:
    input:
        ind_snp_vcf=f"{data_dir}/{{sample}}_snps.vcf"
    output:
        dp_vcf=f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf"
    params:
        maxDP=lambda wildcards: maxDP[samples.index(wildcards.sample)] 
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools --vcf {input.ind_snp_vcf} --maxDP {params.maxDP} --minDP 5 --recode --recode-INFO-all --out {output.dp_vcf}
        mv {output.dp_vcf}.recode.vcf {output.dp_vcf}
        """


# filter for individual depth of invariant files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf", sample=samples)

rule filt_dp_snp_samples:
    input:
        ind_invar_vcf=f"{data_dir}/{{sample}}_invar.vcf"
    output:
        dp_vcf=f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf"
    params:
        maxDP=lambda wildcards: maxDP[samples.index(wildcards.sample)] 
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools --vcf {input.ind_invar_vcf} --maxDP {params.maxDP} --minDP 5 --recode --recode-INFO-all --out {output.dp_vcf}
        mv {output.dp_vcf}.recode.vcf {output.dp_vcf}
        """


samples = ["SRR12424410", "SRR12424411", "SRR12424412", "SRR12424413", "SRR12424416", "SRR12424417", "SRR12424418", "SRR12424419", "SRR12424421", "SRR12424422", "SRR12424423", "SRR3103524"]

# zip snp files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf.gz.tbi", sample=samples)

rule vcf_to_gzvcf_snpfiles:
    input:
        ind_dp_vcf=f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf",
    output:
        ind_dp_gzvcf=f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf.gz",
        ind_dp_gzvcf_tbi=f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf.gz.tbi",
    shell:
        """
        module load  HTSlib/1.18-GCC-12.2.0
        bgzip {input.ind_dp_vcf}
        tabix -p vcf {output.ind_dp_gzvcf}
        cp {output.ind_dp_gzvcf_tbi} {data_dir}/misc
        """

# zip invariant files
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf.gz.tbi", sample=samples)
rule vcf_to_gzvcf_invarfiles:
    input:
        ind_dp_vcf=f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf",
    output:
        ind_dp_gzvcf=f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf.gz",
        ind_dp_gzvcf_tbi=f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf.gz.tbi",
    shell:
        """
        module load HTSlib/1.18-GCC-12.2.0
        bgzip {input.ind_dp_vcf}
        tabix -p vcf {output.ind_dp_gzvcf}
        cp {output.ind_dp_gzvcf_tbi} {data_dir}/misc
        """

rule all:
    input:
        f"{data_dir}/til_caes_snps_filtered_maxdp_mindp5.vcf",
        f"{data_dir}/til_caes_invar_filtered_maxdp_mindp5.vcf"

rule merge_vcfs:
    input:
        snps = expand(f"{data_dir}/{{sample}}_snps_maxdp_mindp5.vcf.gz", sample=samples),
        invar = expand(f"{data_dir}/{{sample}}_invar_maxdp_mindp5.vcf.gz", sample=samples)
    output:
        snps_merged_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp5.vcf",
        invar_merged_vcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp5.vcf"
    shell:
        """
        module load  BCFtools/1.15.1-GCC-11.3.0
        bcftools merge {input.snps}  > {output.snps_merged_vcf}
        bcftools merge {input.invar}  > {output.invar_merged_vcf}
        """

# insert rule that removes hets that do not meet the 25-75 balance ratio

rule filt_hets:

import sys

def update_genotype(gt, ad):
    if gt in ["0/0", "0|0", "1/1", "1|1"] or gt == "./.":
        return gt
    ref_count, alt_count = map(int, ad.split(','))
    total_reads = ref_count + alt_count

    if total_reads == 0:
        return "./."

    allele_ratio = alt_count / total_reads

    if 0.25 <= allele_ratio <= 0.75:
        return gt
    else:
        return "./."

def process_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                # Write header lines to the output file
                outfile.write(line)
            else:
                # Process data lines
                fields = line.strip().split('\t')
                format_field = fields[8].split(':')
                gt_index = format_field.index('GT')
                ad_index = format_field.index('AD')

                # Iterate over samples
                for i in range(9, len(fields)):
                    sample_fields = fields[i].split(':')
                    gt = sample_fields[gt_index]
                    ad = sample_fields[ad_index]

                    # Update genotype
                    updated_gt = update_genotype(gt, ad)
                    sample_fields[gt_index] = updated_gt

                    # Update the line with the modified sample data
                    fields[i] = ':'.join(sample_fields)

                # Write the updated line to the output file
                outfile.write('\t'.join(fields) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input.vcf output.vcf")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    process_vcf(input_file, output_file)
    print(f"Genotype filtering completed. Output written to {output_file}")


# randomly select allele for heterozygous 
import random
import vcf

def randomly_homogenize_genotype(genotype):
    alleles = genotype.split('/')
    if len(alleles) == 2:  # For genotypes in the form 0/1
        random.shuffle(alleles)
        return '/'.join(alleles)
    elif len(alleles) == 2 and '|' in genotype:  # For genotypes in the form 0|1
        random.shuffle(alleles)
        return '|'.join(alleles)
    else:
        return genotype  # For homozygous genotypes, leave unchanged

def modify_vcf(input_vcf, output_vcf):
    vcf_reader = vcf.Reader(open(input_vcf, 'r'))
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)

    for record in vcf_reader:
        for sample in record.samples:
            original_genotype = str(sample['GT'])
            modified_genotype = randomly_homogenize_genotype(original_genotype)
            sample['GT'] = vcf.model._Call(record, sample.sample, modified_genotype)

        vcf_writer.write_record(record)

    vcf_writer.close()

# usage:
input_vcf_file = 'til_caes_snps_filtered_maxdp_mindp5.vcf'
output_vcf_file = 'til_caes_snps_filtered_maxdp_mindp5_het.vcf'
modify_vcf(input_vcf_file, output_vcf_file)


ml BCFtools/1.15.1-GCC-11.3.0
bcftools view -e 'ALT="*"' -O v -o til_caes_snps_filtered_maxdp_mindp5_het_nodel.vcf til_caes_snps_filtered_maxdp_mindp5_het.vcf


#filter for minor allele count and max missing (2/12 samples with a genotype present)
rule mac_filter:
    input:
        filtered_hets_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp5_het_nodel.vcf"
    output:
        filtered_mac_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp5_het_nodel_mm_mac.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools \
            --vcf {input.filtered_hets_vcf} \
            --remove-indels \
            --min-alleles 2 \
            --max-alleles 2 \
            --max-missing-count 10 \
            --mac 3 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mac_vcf}
         mv {output.filtered_mac_vcf}.recode.vcf {output.filtered_mac_vcf}
        """

#filter for max missing (2/12 samples with a genotype present)
rule filter_inv:
    input:
        filtered_invcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp5.vcf"
    output:
        filtered_mm_invcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp5_mm.vcf"
    shell:
        """
        module load VCFtools/0.1.16-GCC-11.2.0
        vcftools \
            --vcf {input.filtered_invcf} \
            --max-missing-count 10 \
            --recode \
            --recode-INFO-all \
            --out {output.filtered_mm_invcf}
         mv {output.filtered_mm_invcf}.recode.vcf {output.filtered_mm_invcf}
        """

#no mac or mm for both inv and snp, test
# bgzip and tabix files
rule vcf_to_gzvcf:
    input:
        var_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp10_hetpy_biallel2.vcf",
        invar_vcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp10.vcf"
    output:
        gz_var_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp10_hetpy_biallel2.vcf.gz",
        tabix_var_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp10_hetpy_biallel2.vcf.gz.tbi",
        gz_invar_vcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp10.vcf.gz",
        tabix_invar_vcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp10.vcf.tbi"
    shell:
        """
        module load HTSlib/1.18-GCC-12.2.0
        bgzip {input.invar_vcf}
        bgzip {input.var_vcf}
        tabix -p vcf {output.gz_invar_vcf}
        tabix -p vcf {output.gz_var_vcf}
        """

rule combine_vcfs:
    input:
       gz_var_vcf=f"{data_dir}/til_caes_snps_filtered_maxdp_mindp5_het_biallel_nodel_macmm.vcf.gz",
       gz_invar_vcf=f"{data_dir}/til_caes_invar_filtered_maxdp_mindp5_mm.vcf.gz"
    output:
       final_vcf=f"{data_dir}/til_caes_allsamples_allsites_final.vcf.gz",
       tabix_final=f"{data_dir}/til_caes_allsamples_allsites_final.vcf.gz.tbi"
    shell:
        """
        module load HTSlib/1.18-GCC-12.2.0
        module load BCFtools/1.15.1-GCC-11.3.0
        bcftools concat {input.gz_var_vcf} {input.gz_invar_vcf} -a -Oz -o {output.final_vcf}
        tabix -p vcf {output.final_vcf}
        """  

