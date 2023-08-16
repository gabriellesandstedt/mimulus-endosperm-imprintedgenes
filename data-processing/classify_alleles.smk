# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm/data/star_pass2"
bed_dir = "/scratch/gds44474/MIMULUS/snps_parents/data"

rule classify_alleles_caes:
    input:
        py_script=f"{star_pass2_dir}/Classify_Alleles.py",
        caes_bed=f"{bed_dir}/final_caes.bed",
        UxT1=f"{star_pass2_dir}/31_S8_STAR_IM62_v3_MD_Split_Q60.bam",
        UxT2=f"{star_pass2_dir}/33_S22_STAR_IM62_v3_MD_Split_Q60.bam",
        UxT3=f"{star_pass2_dir}/48_S28_STAR_IM62_v3_MD_Split_Q60.bam",
        TxU1=f"{star_pass2_dir}/35_S10_STAR_IM62_v3_MD_Split_Q60.bam",
        TxU2=f"{star_pass2_dir}/52_S15_STAR_IM62_v3_MD_Split_Q60.bam",
        TxU3=f"{star_pass2_dir}/45_S14_STAR_IM62_v3_MD_Split_Q60.bam"
    output:
        UxT1_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_1",
        UxT2_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_2",
        UxT3_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_3",
        TxU1_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_1",
        TxU2_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_2",
        TxU3_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_3"
    shell:
        """
        module load Pysam/0.10.0-GCCcore-8.3.0-Python-2.7.16
        python {input.py_script} {input.UxT1} {input.caes_bed} > {output.UxT1_final}
        python {input.py_script} {input.UxT2} {input.caes_bed} > {output.UxT2_final}
        python {input.py_script} {input.UxT3} {input.caes_bed} > {output.UxT3_final}
        python {input.py_script} {input.TxU1} {input.caes_bed} > {output.TxU1_final}
        python {input.py_script} {input.TxU2} {input.caes_bed} > {output.TxU2_final}
        python {input.py_script} {input.TxU3} {input.caes_bed} > {output.TxU3_final}
        """
rule classify_alleles_til:
    input:
        py_script=f"{star_pass2_dir}/Classify_Alleles.py",
        til_bed=f"{bed_dir}/final_til.bed",
        SxL1=f"{star_pass2_dir}/15_S7_STAR_IM62_v3_MD_Split_Q60.bam",
        SxL2=f"{star_pass2_dir}/39_S23_STAR_IM62_v3_MD_Split_Q60.bam",
        SxL3=f"{star_pass2_dir}/46_S26_STAR_IM62_v3_MD_Split_Q60.bam",
        LxS1=f"{star_pass2_dir}/13_S17_STAR_IM62_v3_MD_Split_Q60.bam",
        LxS2=f"{star_pass2_dir}/41_S24_STAR_IM62_v3_MD_Split_Q60.bam",
        LxS3=f"{star_pass2_dir}/50_S30_STAR_IM62_v3_MD_Split_Q60.bam"
    output:
        SxL1_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_1",
        SxL2_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_2",
        SxL3_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_3",
        LxS1_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_1",
        LxS2_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_2",
        LxS3_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_3"
    shell:
        """
        module load Pysam/0.10.0-GCCcore-8.3.0-Python-2.7.16
        python {input.py_script} {input.SxL1} {input.til_bed} > {output.SxL1_final}
        python {input.py_script} {input.SxL2} {input.til_bed} > {output.SxL2_final}
        python {input.py_script} {input.SxL3} {input.til_bed} > {output.SxL3_final}
        python {input.py_script} {input.LxS1} {input.til_bed} > {output.LxS1_final}
        python {input.py_script} {input.LxS2} {input.til_bed} > {output.LxS2_final}
        python {input.py_script} {input.LxS3} {input.til_bed} > {output.LxS3_final}
        """

# assign directories
star_pass2_dir = "/scratch/gds44474/MIMULUS/RNAseq_endosperm_til/data/star_pass2"
bed_dir = "/scratch/gds44474/MIMULUS/snps_parents_til/data"

rule classify_alleles_caes:
    input:
        py_script=f"{star_pass2_dir}/Classify_Alleles.py",
        caes_bed=f"{bed_dir}/final_caes.bed",
        UxT1=f"{star_pass2_dir}/31_S8_STAR_LVR_v1_MD_Split_Q60.bam",
        UxT2=f"{star_pass2_dir}/33_S22_STAR_LVR_v1_MD_Split_Q60.bam",
        UxT3=f"{star_pass2_dir}/48_S28_STAR_LVR_v1_MD_Split_Q60.bam",
        TxU1=f"{star_pass2_dir}/35_S10_STAR_LVR_v1_MD_Split_Q60.bam",
        TxU2=f"{star_pass2_dir}/52_S15_STAR_LVR_v1_MD_Split_Q60.bam",
        TxU3=f"{star_pass2_dir}/45_S14_STAR_LVR_v1_MD_Split_Q60.bam"
    output:
        UxT1_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_1",
        UxT2_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_2",
        UxT3_final=f"{star_pass2_dir}/count_alleles_UTCxTWN_3",
        TxU1_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_1",
        TxU2_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_2",
        TxU3_final=f"{star_pass2_dir}/count_alleles_TWNxUTC_3"
    shell:
        """
        module load Pysam/0.10.0-GCCcore-8.3.0-Python-2.7.16
        python {input.py_script} {input.UxT1} {input.caes_bed} > {output.UxT1_final}
        python {input.py_script} {input.UxT2} {input.caes_bed} > {output.UxT2_final}
        python {input.py_script} {input.UxT3} {input.caes_bed} > {output.UxT3_final}
        python {input.py_script} {input.TxU1} {input.caes_bed} > {output.TxU1_final}
        python {input.py_script} {input.TxU2} {input.caes_bed} > {output.TxU2_final}
        python {input.py_script} {input.TxU3} {input.caes_bed} > {output.TxU3_final}
        """


rule classify_alleles_til:
    input:
        py_script=f"{star_pass2_dir}/Classify_Alleles.py",
        til_bed=f"{bed_dir}/final_til.bed",
        SxL1=f"{star_pass2_dir}/15_S7_STAR_LVR_v1_MD_Split_Q60.bam",
        SxL2=f"{star_pass2_dir}/39_S23_STAR_LVR_v1_MD_Split_Q60.bam",
        SxL3=f"{star_pass2_dir}/46_S26_STAR_LVR_v1_MD_Split_Q60.bam",
        LxS1=f"{star_pass2_dir}/13_S17_STAR_LVR_v1_MD_Split_Q60.bam",
        LxS2=f"{star_pass2_dir}/41_S24_STAR_LVR_v1_MD_Split_Q60.bam",
        LxS3=f"{star_pass2_dir}/50_S30_STAR_LVR_v1_MD_Split_Q60.bam"
    output:
        SxL1_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_1",
        SxL2_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_2",
        SxL3_final=f"{star_pass2_dir}/count_alleles_SOPxLVR_3",
        LxS1_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_1",
        LxS2_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_2",
        LxS3_final=f"{star_pass2_dir}/count_alleles_LVRxSOP_3"
    shell:
        """
        module load Pysam/0.10.0-GCCcore-8.3.0-Python-2.7.16
        python {input.py_script} {input.SxL1} {input.til_bed} > {output.SxL1_final}
        python {input.py_script} {input.SxL2} {input.til_bed} > {output.SxL2_final}
        python {input.py_script} {input.SxL3} {input.til_bed} > {output.SxL3_final}
        python {input.py_script} {input.LxS1} {input.til_bed} > {output.LxS1_final}
        python {input.py_script} {input.LxS2} {input.til_bed} > {output.LxS2_final}
        python {input.py_script} {input.LxS3} {input.til_bed} > {output.LxS3_final}
        """




rule classify_alleles_caes:
    input:
        py_script=f"{star_pass2_dir}/Classify_Alleles.py",
        caes_bed=f"{bed_dir}/final_caes.bed",
        inputs=[
            f"{star_pass2_dir}/31_S8_STAR_IM62_v3_MD_Split_Q60.bam",
            f"{star_pass2_dir}/33_S22_STAR_IM62_v3_MD_Split_Q60.bam",
            f"{star_pass2_dir}/48_S28_STAR_IM62_v3_MD_Split_Q60.bam",
            f"{star_pass2_dir}/35_S10_STAR_IM62_v3_MD_Split_Q60.bam",
            f"{star_pass2_dir}/52_S15_STAR_IM62_v3_MD_Split_Q60.bam",
            f"{star_pass2_dir}/45_S14_STAR_IM62_v3_MD_Split_Q60.bam"
        ]
    output:
        outputs=[
            f"{star_pass2_dir}/count_alleles_UTCxTWN_1.txt",
            f"{star_pass2_dir}/count_alleles_UTCxTWN_2.txt",
            f"{star_pass2_dir}/count_alleles_UTCxTWN_3.txt",
            f"{star_pass2_dir}/count_alleles_TWNxUTC_1.txt",
            f"{star_pass2_dir}/count_alleles_TWNxUTC_2.txt",
            f"{star_pass2_dir}/count_alleles_TWNxUTC_3.txt"
        ]
    run:
        module load Pysam/0.10.0-GCCcore-8.3.0-Python-2.7.16

        for input_file, output_file in zip(input.inputs, output.outputs):
            shell("python {input.py_script} {input_file} {input.caes_bed} > {output_file}")
