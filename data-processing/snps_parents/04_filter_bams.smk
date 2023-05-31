################################################################################
## Add readgroups, mark and remove duplicates, adjust fixmate information, & ensure paired end reads map together
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 04_filter_bams.smk
################################################################################
################################################################################
import os
import pandas as pd
from snakemake.io import expand

# Define the paths to data files
data_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"

# assign sample names of bam files to be filtered
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421"]

# assign all output files to rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_RG.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG.bam.bai", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD.bam.bai", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_NS.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_NS.bam.bai", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam.bai", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_PP.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_PP.bam.bai", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_PP_CS.bam", sample=samples),
        expand(f"{data_dir}/{{sample}}_RG_MD_NS_PP_CS.bam.bai", sample=samples)

# define rule to add or replace read groups
# picard v. 2.27: https://broadinstitute.github.io/picard/
rule add_or_replace_read_groups:
    input:
        bam=f"{data_dir}/{{sample}}_sorted.bam"
    output:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam"
    shell:
        """
        module load picard/2.27.4-Java-13.0.2
        echo -e "\\n["$(date)"]\\n Add read groups..\\n"
        java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
            I={input.bam} \
            O={output.RG_bam} \
            RGID={wildcards.sample} \
            RGLB=lib_{wildcards.sample} \
            RGPL=illumina \
            RGPU=unit_{wildcards.sample} \
            RGSM={wildcards.sample}
        """

# define rule to index bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule index_RG_bam:
    input:
        RG_bam=f"{data_dir}/{{sample}}_RG.bam"
    output:
        RG_bai=f"{data_dir}/{{sample}}_RG.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.RG_bam} {output.RG_bai}
        """

# define rule to mark and remove duplicates 
# picard v. 2.27: https://broadinstitute.github.io/picard/
rule mark_duplicates:
    input:
        RG_bam=f"{data_dir}/{{sample}}_sorted_RG.bam"
    output:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam",
        MD_log=f"{data_dir}/{{sample}}_RG_MDlog.txt"
    shell:
        """
        module load picard/2.27.4-Java-13.0.2
        echo -e "\\n["$(date)"]\\n mark and remove duplicates..\\n"
        java -jar $EBROOTPICARD/picard.jar MarkDuplicates I={input.RG_bam} O={output.MD_bam} REMOVE_DUPLICATES=TRUE M={output.MD_log}
        """

# define rule to index bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule index_MD_bam:
    input:
        MD_bam=f"{data_dir}/{{sample}}_RG_MD.bam"
    output:
        MD_bai=f"{data_dir}/{{sample}}_RG_MD.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.MD_bam} {output.MD_bai}
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

# define rule to index bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule index_NS_bam:
    input:
        NS_bam=f"{data_dir}/{{sample}}_RG_MD_NS.bam"
    output:
        NS_bai=f"{data_dir}/{{sample}}_RG_MD_NS.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.NS_bam} {output.NS_bai}
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

# define rule to index bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule index_FM_bam:
    input:
        FM_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam"
    output:
        FM_bai=f"{data_dir}/{{sample}}_RG_MD_NS_FM.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.FM_bam} {output.FM_bai}
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

# define rule to index bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule index_PP_bam:
    input:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    output:
        PP_bai=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.PP_bam} {output.PP_bai}
        """       

# define rule to sort bam files by coordinates
# samtools v 1.16: https://github.com/samtools/samtools
rule coord_sort:
    input:
        PP_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP.bam"
    output:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_NS_FM_PP_CS.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n sort bam by coordinate...\\n"
        samtools sort -o {output.CS_bam} {input.PP_bam}
        """

# define rule to index bam file
# samtools v 1.16: https://github.com/samtools/samtools
rule index_CS_bam:
    input:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_FM_PP_CS.bam"
    output:
        CS_bai=f"{data_dir}/{{sample}}_RG_MD_FM_PP_CS.bam.bai"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\n Index BAM file...\\n"
        samtools index {input.CS_bam} {output.CS_bai}
        """

# define rule that details bam files
# samtools v 1.16: https://github.com/samtools/samtools
rule assess_quality:
    input:
        CS_bam=f"{data_dir}/{{sample}}_RG_MD_FM_PP_CS.bam"
    shell:
        """
        module load SAMtools/1.16.1-GCC-11.3.0
        echo -e "\\n["$(date)"]\\ run samtools flagstat on final bam file...\\n"
        samtools flagstat {input.CS_bam}
        """
