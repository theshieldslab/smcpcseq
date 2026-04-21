#!/bin/bash

#SBATCH --array=1-40
#SBATCH --job-name=fastqc
#SBATCH --partition=hpg-default
#SBATCH --account=robert.shields
#SBATCH --qos=robert.shields-b
#SBATCH --mail-user=ariverosw@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --output stdout_%A_%a_fastqc.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=8gb
#SBATCH --time=03:00:00

start=`date +%s`

### Changing directory
cd /blue/robert.shields/arwalker/fastqc

### SampleID
sampleID=$( eval cat ../samples.txt | sed -n "${SLURM_ARRAY_TASK_ID},${SLURM_ARRAY_TASK_ID}p" | cut -f 1)
#sampleID=${SLURM_ARRAY_TASK_ID}

### Running FASTQC
ml fastqc
fastqc -t ${SLURM_CPUS_PER_TASK} /blue/robert.shields/anasolanomorales/tn_seq/samples/${sampleID}.fastq.gz --outdir ./

### Print time
end=`date +%s`
runtime=$((end-start))
out=$(echo "scale=2; (${runtime}/1)" | bc); echo "${out} seconds"
