################################################################################
## create bed files for tilingii and caespitosa with 6 columns:
## 1|CHROMOSOME 2|POSITION START 3|POSITION END 4|SNP INDIVIDUAL 1 5|SNP INDIVIDUAL 2 6|GENE NAME
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete --latency-wait 60 --cores 4 -s 07_process_parent_bed_files.smk
################################################################################
################################################################################

# assign directories
data_dir = "/scratch/gds44474/MIMULUS/snps_parents_til/data"
ref_dir = "/scratch/gds44474/MIMULUS/ref_genome_til"

# assign files 
gff3= "MtilingiivarLVRv1.1.primaryTrs.gff3"


# assign rule for all expected output files from this script 
rule all:
    input:
        f"{data_dir}/genes.bed",
        f"{data_dir}/final_caes.bed",
        f"{data_dir}/final_til.bed",
        f"caes_genes_no.txt",
        f"til_genes_no.txt"
        
# assign rule to convert gff3 file to bed file 
rule gff_to_bed:
    input:
        gff=f"{ref_dir}/{gff3}"
    output:
        genes_bed=f"{data_dir}/genes.bed"
    shell:
        """
        module load BEDOPS/2.4.39-foss-2019b
        gff2bed < {input.gff} > {output.genes_bed}
        """ 

# UTC1 SNP is in column 4
# TWN36 SNP is in column 5
rule modify_caes_bed:
    input:
        genes_bed=f"{data_dir}/genes.bed",
        caes_bed=f"{data_dir}/caes.bed"
    output:
        final_caes_bed=f"{data_dir}/final_caes.bed"
    shell:
        """
        module load BEDTools/2.30.0-GCC-10.2.0

        awk 'BEGIN{{OFS="\t"}} {{
            if (FNR == 0) print
            else {{
                for (i=11; i<=12; i++) {{
                    split($i, a, ":")
                    $i = a[1]
                    if ($i == "0/0" || $i == "0|0") $i = $6
                    if ($i == "1/1" || $i == "1|1") $i = $7
                }}
                print
            }}
        }}' {input.caes_bed} |
        bedtools intersect -a stdin -b {input.genes_bed} -wb |
        awk -v OFS="\t" '$20 ~ /gene/ {{ print $1,$2,$3,$11,$12,$16 }}' > {output.final_caes_bed}
        """

# SOP12 SNP is in column 4
# LVR1 SNP is in column 5
rule modify_til_bed:
    input:
        genes_bed=f"{data_dir}/genes.bed",
        til_bed=f"{data_dir}/til.bed"
    output:
        final_til_bed=f"{data_dir}/final_til.bed"
    shell:
        """
        module load BEDTools/2.30.0-GCC-10.2.0

        awk 'BEGIN{{OFS="\t"}} {{
            if (FNR == 0) print
            else {{
                for (i=11; i<=12; i++) {{
                    split($i, a, ":")
                    $i = a[1]
                    if ($i == "0/0" || $i == "0|0") $i = $6
                    if ($i == "1/1" || $i == "1|1") $i = $7
                }}
                print
            }}
        }}' {input.til_bed} |
        bedtools intersect -a stdin -b {input.genes_bed} -wb |
        awk -v OFS="\t" '$20 ~ /gene/ {{ print $1,$2,$3,$11,$12,$16 }}' > {output.final_til_bed}
        """

# caespitosa: count number of genes with a snp
rule count_snp_caes:
    input:
        final_caes_bed=f"{data_dir}/final_caes.bed"
    output:
        genes_no=f"{data_dir}/caes_genes_no.txt"
    shell:
        """
        cut -f 6  {input.final_caes_bed} | sort | uniq | wc -l > {output.genes_no}
        """ 

# tilingii: count number of genes with a snp
rule count_snp_til:
    input:
        final_til_bed=f"{data_dir}/final_til.bed"
    output:
        genes_no=f"{data_dir}/til_genes_no.txt"
    shell:
        """
        cut -f 6  {input.final_til_bed} | sort | uniq | wc -l > {output.genes_no}
        """ 
