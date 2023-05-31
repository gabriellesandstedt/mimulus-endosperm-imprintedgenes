################################################################################
## Mask repetitive elements of reference genome, index and align with STAR
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 02_STAR_align.smk
################################################################################
################################################################################
import os

# assign directories
ref_dir="/scratch/gds44474/MIMULUS/ref_genome"
repeatmasker_dir="/scratch/gds44474/MIMULUS/ref_genome/RepeatMasker"
STAR_genome_dir="/scratch/gds44474/MIMULUS/ref_genome/RepeatMasker/STAR_genome"
data_dir="/scratch/gds44474/MIMULUS/RNAseq_endosperm/data"
star_pass2_dir="/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"

# assign genome files
ref = "Mimulus_guttatus_var_IM62.mainGenome.fasta"
masked_ref = "Mimulus_guttatus_var_IM62.mainGenome.fasta.masked"
masked_ref2 = "Mimulus_guttatus_var_IM62.mainGenome.masked.fasta"
gff = "MguttatusvarIM62v3.1.primaryTrs.gff3"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]

# define all output files to rule all
rule all:
    input:
        f"{repeatmasker_dir}/{masked_ref2}",
        expand(f"{data_dir}/{{sample}}.bam", sample=samples)

# define rule to mask repetitive elements in the genome
rule repeatmasker:
    input:
        ref_genome = f"{ref_dir}/{ref}"
    output:
        repmask_dir = directory(repeatmasker_dir)
    shell:
        """
        ml RepeatMasker/4.1.2-p1-foss-2020b
        RepeatMasker {input.ref_genome} -species Erythranthe -dir {output.repmask_dir}
        """

# define rule to rename the masked genome
rule rename_masked_fasta:
    input:
        masked_fa = f"{repeatmasker_dir}/{masked_ref}"
    output:
        masked_fa2 = f"{repeatmasker_dir}/{masked_ref2}"
    shell:
        """
        cp {input.masked_fa} {output.masked_fa2}
        """

# define rule to index the masked genome for STAR alignment
rule create_star_index:
    input:
        masked_fa2 = f"{repeatmasker_dir}/{masked_ref2}",
        gff = f"{ref_dir}/{gff}"
    output:
        genome_dir = directory(STAR_genome_dir)
    shell:
        """
        ml STAR/2.7.10a-GCC-8.3.0
        STAR --runThreadN 12 --runMode genomeGenerate --genomeDir {output.genome_dir} --genomeFastaFiles {input.masked_fa2} --sjdbGTFfile {input.gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --genomeSAindexNbases 13 --sjdbGTFtagExonParentGene Parent
        """

# define rule to align fastqs to masked reference genome
rule star_alignment_pass1:
    input:
        masked_fa2 = f"{repeatmasker_dir}/{masked_ref2}",
        gff = f"{ref_dir}/{gff}",
        R1 = expand(f"{data_dir}/{{sample}}_R1.trim.fastq.gz", sample=samples),
        R2 = expand(f"{data_dir}/{{sample}}_R2.trim.fastq.gz", sample=samples),
        genome_dir = directory(STAR_genome_dir)
    params:
        sample = "{sample}"
    output:
        bam = f"{data_dir}/{{sample}}.bam"
    shell:
        """
        ml STAR/2.7.10a-GCC-8.3.0
        STAR --runThreadN 12 --genomeDir {input.genome_dir} --sjdbGTFfile {input.gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --alignIntronMin 20 --alignIntronMax 10000 --outFilterMismatchNoverReadLmax 0.05 --outSAMmapqUnique 60 --outFileNamePrefix {output.bam} --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMattrRGline ID:{params.sample} LB:IM62.v3 DS:RNAseq PU:NovaSeq6000 PL:Illumina SM:{params.sample}
        """

# define rule to find SJ.out.tab files
rule find_sj_files:
    input:
        data_dir = data_dir
    output:
        sj_files = f"{data_dir}/sj_files.txt"
    shell:
        """
        find {input.data_dir} -name "*SJ.out.tab" | tr '\n' ' ' > {output.sj_files}
        """

# define rule to align fastqs to masked reference genome using SJ.out.tab files, which contain info on split junctions 
# this is a second pass at the STAR alignment to improve alignment quality 
rule star_alignment_pass2:
    input:
        masked_fa2 = f"{repeatmasker_dir}/{masked_ref2}",
        gff = f"{ref_dir}/{gff}",
        R1 = expand(f"{data_dir}/{{sample}}_R1.trim.fastq.gz", sample=samples),
        R2 = expand(f"{data_dir}/{{sample}}_R2.trim.fastq.gz", sample=samples),
        sj_files = output.sj_files,
        genome_dir = directory(STAR_genome_dir)
    params:
        sample = "{sample}"
    output:
        directory = directory(star_pass2_dir),
        bam = f"{data_dir}/{{sample}}.bam"
    shell:
        """
        ml STAR/2.7.10a-GCC-8.3.0
        STAR --runThreadN 12 --genomeDir {input.genome_dir} --sjdbFileChrStartEnd $(cat {input.sj_files}) --sjdbGTFfile {input.gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} --alignIntronMin 20 --alignIntronMax 10000 --outFilterMismatchNoverReadLmax 0.05 --outFileNamePrefix {output.bam} --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outSAMattributes All --outSAMattrRGline ID:{params.sample} LB:IM62.v3 DS:RNAseq
        """



