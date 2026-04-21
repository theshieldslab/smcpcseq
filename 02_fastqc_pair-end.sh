#!/usr/bin/bash

#SBATCH --array=1-40
#SBATCH --job-name=bar#split
#SBATCH --partition=hpg-default
#SBATCH --account=robert.shields
#SBATCH --mail-user=ariverosw@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --output stdout_%A_%a_barcodesplit.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=04:00:00

### Timestamp
start=`date +%s`

### Changing directory
cd /blue/robert.shields/arwalker/samples

### Variable creation
sampleID=$( cat ../barcodes.txt | sed -n "${SLURM_ARRAY_TASK_ID},${SLURM_ARRAY_TASK_ID}p" | cut -f 1 )
barcode=$( cat ../barcodes.txt | sed -n "${SLURM_ARRAY_TASK_ID},${SLURM_ARRAY_TASK_ID}p" | cut -f 2 )

echo -e ${SLURM_ARRAY_TASK_ID}"\t"${sampleID}"\t"${barcode}

### Samples splitting
echo "Splitting fastq file"
zcat /blue/robert.shields/anasolanomorales/fastq/Tn_1_S0_R1_001.fastq.gz | tr ' ' '_' | paste - - - - | awk -v sample=${sampleID} -v code=${barcode} '{ if( substr($2,1,8)==code ) sum=sum+1 ; if( substr($2,1,8)==code ) print "@"sample"_"code"_"sum"\n"substr($2,9,16)"\n"$3"\n"substr($4,9,16) }' > ${sampleID}.fastq

### Compressing fastq2gz
echo "fastq2gz"
gzip ${sampleID}.fastq

### Counting reads
count=$( zcat ${sampleID}.fastq.gz | tail -4 | head -1 | tr '_' '\t' | cut -f 3 )
echo -e "Reads count = "${count}

### Print time
end=`date +%s`
runtime=$((end-start))
out=$(echo "scale=2; (${runtime}/1)" | bc); echo "${out} seconds"