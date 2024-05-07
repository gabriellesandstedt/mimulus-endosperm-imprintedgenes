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
