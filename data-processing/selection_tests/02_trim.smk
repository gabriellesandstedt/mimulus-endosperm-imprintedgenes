################################################################################
## check quality of fastqs and trim adapter sequences
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: snakemake/7.22.0-foss-2022a
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 13 -s 02_trim.smk
################################################################################
################################################################################
import os

# assign data directory
data_dir = "/scratch/gds44474/MIMULUS/selection_tests/data"
qc1_dir = f"{data_dir}/FASTQC1"
qc2_dir = f"{data_dir}/FASTQC2"
log_dir = f"{data_dir}/logs"

# assign samples to be trimmed
samples = ["SRR12424410", "SRR3103524", "SRR12424419", "SRR12424421", "SRR12424411", "SRR12424412", "SRR12424417", "SRR12424423", "SRR12424422", "SRR12424418","SRR12424416", "SRR12424413"]

# define output files for rule all
rule all:
    output:
        expand(f"{qc1_dir}/{{sample}}_1_fastqc.html", sample=samples),
        expand(f"{qc1_dir}/{{sample}}_2_fastqc.html", sample=samples),
        expand(f"{data_dir}/{{sample}}_1_trim.fq.gz", sample=samples),
        expand(f"{data_dir}/{{sample}}_2_trim.fq.gz", sample=samples),
        expand(f"{qc2_dir}/{{sample}}_1_trim_fastqc.html", sample=samples),
        expand(f"{qc2_dir}/{{sample}}_2_trim_fastqc.html", sample=samples)


# define rule to assess quality of fastqs with FastQC
# FASTQC v 0.11.9 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_raw:
    input:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    output:
        fqc1=f"{qc1_dir}/{{sample}}_1_fastqc.html",
        fqc2=f"{qc1_dir}/{{sample}}_2_fastqc.html"
    log:
        log1=f"{log_dir}/{{sample}}_1.log",
        log2=f"{log_dir}/{{sample}}_2.log"
    shell:
        """
        module load FastQC/0.11.9-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq1} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq1} &> {log.log1} 
        echo -e "\\n["$(date)"]\\n FastQC round 1, read 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.fq2} ...\\n"
        fastqc -o {qc1_dir} --noextract {input.fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 1, read 2 finished ...\\n"
        """

# define rule to trim adapter sequences and filter raw fastq reads using trimmomatic
# trimmomatic v 0.39: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
rule trimmomatic:
    input:
        fq1=f"{data_dir}/{{sample}}_1.fastq.gz",
        fq2=f"{data_dir}/{{sample}}_2.fastq.gz"
    output:
        trim_fq1=f"{data_dir}/{{sample}}_1_trim.fq.gz",
        trim_fq2=f"{data_dir}/{{sample}}_1_trim_unpaired.fq.gz",
        trim_fq3=f"{data_dir}/{{sample}}_2_trim.fq.gz",
        trim_fq4=f"{data_dir}/{{sample}}_2_trim_unpaired.fq.gz"
    shell:
        """
        module load Trimmomatic/0.39-Java-13
        echo -e "\\n["$(date)"]\\n Run Trimmomatic on raw fastq file {input.fq1} and {input.fq2} ...\\n"
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 4 -phred33 {input.fq1} {input.fq2} {output.trim_fq1} {output.trim_fq2} {output.trim_fq3} {output.trim_fq4} ILLUMINACLIP:TruSeq3.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
        echo -e "\\n["$(date)"]\\n Trimmomatic finished ...\\n"
        """
# define rule to assess quality of trimmed fastqs with FastQC
# FASTQC v 0.11.9 : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
rule fastqc_trimmed:
    input:
        trim_fq1=f"{data_dir}/{{sample}}_1_trim.fq.gz",
        trim_fq2=f"{data_dir}/{{sample}}_2_trim.fq.gz"
    output:
        trim_fqc1=f"{qc2_dir}/{{sample}}_1_trim_fastqc.html",
        trim_fqc2=f"{qc2_dir}/{{sample}}_2_trim_fastqc.html"
    log:
        log1=f"{log_dir}/trim_{{sample}}_1.log",
        log2=f"{log_dir}/trim_{{sample}}_2.log"
    shell:
        """
        module load FastQC/0.11.9-Java-11
        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq1} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.trim_fq1} &> {log.log1}
        echo -e "\\n["$(date)"]\\n FastQC round 2, read 1 finished ...\\n"

        echo -e "\\n["$(date)"]\\n Run FastQC on fastq file {input.trim_fq2} ...\\n"
        fastqc -o {qc2_dir} --noextract {input.trim_fq2} &> {log.log2}
        echo -e "\\n["$(date)"]\\n FastQC round 2, read 2 finished ...\\n"
        """
