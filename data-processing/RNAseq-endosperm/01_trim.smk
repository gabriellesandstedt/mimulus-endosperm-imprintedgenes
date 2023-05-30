################################################################################
## concatenate raw fastqs across lanes and trim fastqs
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 01_trim.smk
################################################################################
################################################################################
import os

# assign directories
raw_dir = "/scratch/gds44474/MIMULUS/rna_seq/rawfiles"
data_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm/data"
qc1_dir = f"{data_dir}/FASTQC1"
qc2_dir = f"{data_dir}/FASTQC2"
log_dir = f"{data_dir}/logs"

samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11","34_S9","53_S16","49_S29"]

# define output files for rule all
rule all:
    input:
        expand(f"{data_dir}/{{sample}}_R1.fastq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_R1.fastq.gz", sample=samples),
        expand(f"{qc1_dir}/{{sample}}_R1_fastqc.html", sample=samples),
        expand(f"{qc1_dir}/{{sample}}_R2_fastqc.html", sample=samples),
        expand(f"{data_dir}/{{sample}}_R1_trim.fq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_R2_trim.fq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_R1_trim_fastqc.html", sample=samples),
        expand(f"{data_dir}/{{sample}}_R2_trim_fastqc.html", sample=samples)

# define rule to combine raw fastqs across lanes 
rule combine_fastq:
    input:
        r1L1=f"{raw_dir}/{{sample}}_L001_R1_001.fastq.gz",
        r1L2=f"{raw_dir}/{{sample}}_L002_R1_001.fastq.gz",
        r2L1=f"{raw_dir}/{{sample}}_L001_R2_001.fastq.gz",
        r2L2=f"{raw_dir}/{{sample}}_L002_R2_001.fastq.gz"
    output:
        R1=f"{data_dir}/{{sample}}_R1.fastq.gz",
        R2=f"{data_dir}/{{sample}}_R2.fastq.gz"
    shell:
        """
        gunzip -c {input.r1L1} {input.r1L2} | gzip > {output.R1}
        gunzip -c {input.r2L1} {input.r2L2} | gzip > {output.R2}
        """

# define rule to assess quality of fastqs with FastQC
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_raw:
    input:
        fq1=f"{data_dir}/{{sample}}_R1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_R2.fastq.gz"
    output:
        fqc1=f"{qc1_dir}/{{sample}}_R1_fastqc.html",
        fqc2=f"{qc1_dir}/{{sample}}_R2_fastqc.html"
    log:
        log1=f"{log_dir}/fqc_{{sample}}_R1.log",
        log2=f"{log_dir}/fqc_{{sample}}_R2.log"
    shell:
        """
        module load FastQC/0.12.1-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq1} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq1} &> {log.log1}
        echo -e "\\n["$(date)"]\\n FastQC round 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq2} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 1 finished ...\\n"
        """

# define rule to trim adapter sequences and filter raw fastq reads using trimmomatic
# trimmomatic v 0.39: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
rule trimmomatic:
    input:
        fq1=f"{data_dir}/{{sample}}_R1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_R2.fastq.gz"
    output:
        trim_fq1=f"{data_dir}/{{sample}}_R1_trim.fq.gz",
        trim_fq2=f"{data_dir}/{{sample}}_R1_trim_unpaired.fq.gz",
        trim_fq3=f"{data_dir}/{{sample}}_R2_trim.fq.gz",
        trim_fq4=f"{data_dir}/{{sample}}_R2_trim_unpaired.fq.gz"
    shell:
        """
        module load Trimmomatic/0.39-Java-1.8.0_144
        echo -e "\\n["$(date)"]\\n Run Trimmomatic on raw fastq file {input.fq1} and {input.fq2} ...\\n"
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 -phred33 {input.fq1} {input.fq2} {output.trim_fq1} {output.trim_fq2} {output.trim_fq3} {output.trim_fq4} ILLUMINACLIP:TruSeq3.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
        echo -e "\\n["$(date)"]\\n Trimmomatic finished ...\\n"
        """
# define rule to assess quality of trimmed fastqs with FastQC
# FASTQC v 0.12.1 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_trimmed:
    input:
        trim_fq1=f"{data_dir}/{{sample}}_R1.fastq.gz",
        trim_fq2=f"{data_dir}/{{sample}}_R2.fastq.gz"
    output:
        trim_fqc1=f"{qc2_dir}/{{sample}}_R1_fastqc.html",
        trim_fqc2=f"{qc2_dir}/{{sample}}_R2_fastqc.html"
    log:
        log1=f"{log_dir}/trim_fqc_{{sample}}_R1.log",
        log2=f"{log_dir}/trim_fqc_{{sample}}_R2.log"
    shell:
        """
        module load FastQC/0.12.1-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq1} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.trim_fq1} &> {log.log1}
        echo -e "\\n["$(date)"]\\n FastQC round 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq2} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.trim_fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 1 finished ...\\n"
        """
