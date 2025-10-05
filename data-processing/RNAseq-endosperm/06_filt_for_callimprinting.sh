#!/bin/bash
#SBATCH --partition=batch
#SBATCH --time=100:00:00
#SBATCH --nodes=2 #Note: HPCC nodes have 20+ CPU cores each
#SBATCH --ntasks=4 #parallel tasks/processes run on separate CPU cores or nodes (coarse-grained internode/intercore parallelization via MPI; requires OpenMPI/impi) Note: sbatch does not launch the tasks, only requests resources and submits the batch script. ntasks tells Slurm that N parallel tasks will be launched and to allocate resources accordingly. Parallel tasks are launched by the script.
#SBATCH --mem=72G
#SBATCH --job-name cigar_filt
#SBATCH --output 06.out
#SBATCH --error 06.error


files=(
  UTC1xTWN36_3_STAR_LVR_v1_MD_Split_Q60.bam
  UTC1xTWN36_2_STAR_LVR_v1_MD_Split_Q60.bam
  UTC1xTWN36_1_STAR_LVR_v1_MD_Split_Q60.bam
  TWN36xUTC1_3_STAR_LVR_v1_MD_Split_Q60.bam
  TWN36xUTC1_2_STAR_LVR_v1_MD_Split_Q60.bam
  TWN36xUTC1_1_STAR_LVR_v1_MD_Split_Q60.bam
  SOP12xLVR1_3_STAR_LVR_v1_MD_Split_Q60.bam
  SOP12xLVR1_2_STAR_LVR_v1_MD_Split_Q60.bam
  SOP12xLVR1_1_STAR_LVR_v1_MD_Split_Q60.bam
  LVR1xSOP12_3_STAR_LVR_v1_MD_Split_Q60.bam
  LVR1xSOP12_2_STAR_LVR_v1_MD_Split_Q60.bam
  LVR1xSOP12_1_STAR_LVR_v1_MD_Split_Q60.bam
)

for bam in "${files[@]}"; do
  out="${bam%.bam}.noSH.bam"
  echo "[*] Filtering $bam -> $out"
  samtools view -h -@ "$threads" "$bam" \
    | awk 'BEGIN{OFS="\t"} /^@/ || $6 !~ /[SH]/' \
    | samtools view -@ "$threads" -b -o "$out" -
  samtools index -@ "$threads" "$out"
done

echo "Done."
