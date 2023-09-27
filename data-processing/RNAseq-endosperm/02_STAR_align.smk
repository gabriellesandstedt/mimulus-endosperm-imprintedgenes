################################################################################
## Mask repetitive elements of reference genome, index and align with STAR
################################################################################
################################################################################
## AUTHOR: Gabrielle Sandstedt
## snakemake version: 6.3.0
## command to run snakemake script: snakemake --rerun-incomplete  --latency-wait 60  --cores 4 -s 02_STAR_align.smk
################################################################################
################################################################################
import os

# assign directories
ref_dir="/scratch/gds44474/MIMULUS/ref_genome_til"
repeatmasker_dir="/scratch/gds44474/MIMULUS/ref_genome_til/RepeatMasker"
STAR_genome_dir="/scratch/gds44474/MIMULUS/ref_genom_tile/RepeatMasker/STARgenome"
data_dir="/scratch/gds44474/MIMULUS/RNAseq_endosperm/data"
star_pass2_dir="/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"

# assign genome files
ref = "Mimulus_tilingii_var_LVR.mainGenome.fasta"
masked_ref = "Mimulus_tilingii_var_LVR.mainGenome.fasta.masked"
masked_ref2 = "Mimulus_tilingii_var_LVR.mainGenome.masked.fasta"
gff = "MtilingiivarLVRv1.1.primaryTrs.gff3"

# assign samples
samples = ["13_S17", "41_S24", "50_S30", "15_S7", "39_S23", "46_S26", "35_S10", "52_S15", "45_S14", "31_S8", "33_S22", "48_S28", "44_S13", "47_S27", "32_S21", "36_S11", "34_S9", "53_S16", "49_S29"]

# define all output files to rule all
rule all:
    input:
        f"{repeatmasker_dir}/{masked_ref2}",
        expand(f"{data_dir}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=samples),
        expand(f"{star_pass2_dir}/{{sample}}.Aligned.sortedByCoord.out.bam", sample=samples)

# define rule to mask repetitive elements in the genome
# repeat masker v 4.1.2 : https://www.repeatmasker.org
# I first ran this rule separately
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

# define rule to rename the masked genome ; this help with downstream analyses
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
# STAR v 2.7 : https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
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
# STAR v 2.7 : https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
rule star_alignment_pass1:
    input:
        gff=f"{ref_dir}/{gff}",
        rd1=f"{data_dir}/{{sample}}_R1_trim.fastq.gz",
        rd2=f"{data_dir}/{{sample}}_R2_trim.fastq.gz"
    output:
        bam=f"{data_dir}/{{sample}}.Aligned.sortedByCoord.out.bam",
        log_final=f"{data_dir}/{{sample}}.Log.final.out",
        sj=f"{data_dir}/{{sample}}.SJ.out.tab"
    shell:
        """
        ml STAR/2.7.10a-GCC-8.3.0
        STAR --runThreadN 12 --genomeDir {STAR_genome_dir} --sjdbGTFfile {input.gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --readFilesCommand zcat --readFilesIn {input.rd1} {input.rd2} --alignIntronMin 20 --alignIntronMax 10000 --outFilterMismatchNoverReadLmax 0.05 --outSAMmapqUnique 60 --outFileNamePrefix {data_dir}/{wildcards.sample}. --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMattrRGline ID:{wildcards.sample} LB:IM62.v3 DS:RNAseq PU:NovaSeq6000 PL:Illumina SM:{wildcards.sample}
        """

# define rule to align fastqs to masked reference genome using SJ.out.tab files, which contain info on split junctions 
# this is a second pass at the STAR alignment to improve alignment quality 
# STAR v 2.7 : https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
rule star_alignment_pass2:
    input:
        gff = f"{ref_dir}/{gff}",
        rd1 = f"{data_dir}/{{sample}}_R1_trim.fastq.gz",
        rd2 = f"{data_dir}/{{sample}}_R2_trim.fastq.gz",
        sj_files = f"{data_dir}/{{sample}}.SJ.out.tab"
    output:
        bam2=f"{star_pass2_dir}/{{sample}}.Aligned.sortedByCoord.out.bam",
        log_final2=f"{star_pass2_dir}/{{sample}}.Log.final.out",
        sj2=f"{star_pass2_dir}/{{sample}}.SJ.out.tab"
    shell:
        """
        ml STAR/2.7.10a-GCC-8.3.0
        STAR --runThreadN 12 --genomeDir {STAR_genome_dir} --sjdbFileChrStartEnd {input.sj_files} --sjdbGTFfile {input.gff} --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --readFilesCommand zcat --readFilesIn {input.rd1} {input.rd2} --alignIntronMin 20 --alignIntronMax 10000 --outFilterMismatchNoverReadLmax 0.05 --outFileNamePrefix {star_pass2_dir}/{wildcards.sample}. --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 --outSAMattributes All --outSAMattrRGline ID:{wildcards.sample} LB:IM62.v3 DS:RNAseq
        """

