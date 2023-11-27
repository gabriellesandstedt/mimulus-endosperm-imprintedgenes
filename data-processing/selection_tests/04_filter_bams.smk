################################################################################
## Add readgroups, mark and remove duplicates, adjust fixmate information, & ensure paired end reads map together
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: snakemake/7.22.0-foss-2022a
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60 -s 04_filter_bams.smk
################################################################################
################################################################################
import os
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/snps_parents_til/data"

# assign sample names of bam files to be filtered
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421","SRR12424411", "SRR12424412","SRR12424417","SRR12424423", "SRR12424422", "SRR12424418","SRR12424416", "SRR12424413"]

# assign all output files to rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_RG.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_bamQC.pdf", sample=samples)

# define rule to add or replace read groups
# picard v. 2.27: https://broadinstitute.github.io/picard/
# samtools v 1.16: https://github.com/samtools/samtools
rule add_or_replace_read_groups:
    input:
        bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam",
        RG_bai=f"{data_dir}/{{sample}}_RG.bam.bai"
    shell:
        """
        module load picard/2.27.5-Java-15
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Add read groups..\\n"
        java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.RG_bam} \
            RGID={wildcards.sample} \
            RGLB=lib_{wildcards.sample} \
            RGPL=illumina \
            RGPU=unit_{wildcards.sample} \
            RGSM={wildcards.sample}   
        samtools index {output.RG_bam}
        """

# define rule to mark and remove duplicates 
# picard v. 2.27: https://broadinstitute.github.io/picard/
# samtools v 1.16: https://github.com/samtools/samtools
rule mark_duplicates:
    input:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam"
    output:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam",
        MD_log=f"{data_dir}/{{sample}}_RG_MDlog.txt",
        MD_bai=f"{data_dir}/{{sample}}_RG_MD.bam.bai"
    shell:
        """
        module load picard/2.27.5-Java-15
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n mark and remove duplicates..\\n"
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={input.RG_bam} O={output.MD_bam} REMOVE_DUPLICATES=TRUE M={output.MD_log}
        samtools index {output.MD_bam}
        """

# define rule to sort bam files by name
# samtools v 1.16: https://github.com/samtools/samtools
rule namesort:
    input:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam"
    output:
        NS_bam=f"{data_dir}/{{sample}}_RG_MD_NS.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n sort bam by name...\\n"
        samtools sort -o {output.NS_bam} -n {input.MD_bam}
        """
        
# define rule to correct or adjust fixmate information of bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule fixmate:
    input:
        NS_bam=f"{data_dir}/{{sample}}_RG_MD_NS.bam"
    output:
        FM_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n run samtools fixmate...\\n"
        samtools fixmate -r {input.NS_bam} {output.FM_bam}
        """
        
# define rule to ensure that paired end reads map together
# samtools v 1.16: https://github.com/samtools/samtools
rule proper_pair:
    input:
        FM_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam"
    output:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n ensure proper pairing...\\n"
        samtools view -b -f 2 -F 2048 {input.FM_bam} > {output.PP_bam}
        """    

# define rule to sort bam files by coordinates
# samtools v 1.16: https://github.com/samtools/samtools
rule coord_sort:
    input:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    output:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam",
        CS_bai=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n sort bam by coordinate...\\n"
        samtools sort -o {output.CS_bam} {input.PP_bam}
        samtools index {output.CS_bam}
        """
        
# define rule to determine the quality control of alignment sequencing data
# qualimap v 2.2.1: http://qualimap.conesalab.org
rule qualimap:
    input:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam"
    output:
        quali_bam=f"{data_dir}/{{sample}}_bamQC.pdf"
    shell:
        """
        module load Qualimap/2.2.1-foss-2021b-R-4.1.2
        qualimap bamqc -bam {input.CS_bam} -outfile {output.quali_bam} 
        """

