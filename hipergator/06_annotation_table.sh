
#!/usr/bin/bash

#SBATCH --array=1-40
#SBATCH --job-name=mapping
#SBATCH --partition=hpg-default
#SBATCH --account=robert.shields
#SBATCH --mail-user=ariverosw@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --output stdout_%A_%a_mapping.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=04:00:00

### Timestamp
start=`date +%s`

### Changing directory
cd /blue/robert.shields/arwalker/mapping

### Variable creation
sampleID=$( cat ../samples.txt | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f 1 )

echo -e ${SLURM_ARRAY_TASK_ID}"\t"${sampleID}

### Mapping
echo "Mapping fastq file"
ml bwa

bwa aln /blue/robert.shields/anasolanomorales/references/UA159.fasta ../samples/${sampleID}.fastq.gz > ${sampleID}_raw_aln.sai
bwa samse /blue/robert.shields/anasolanomorales/references/UA159.fasta ${sampleID}_raw_aln.sai ../samples/${sampleID}.fastq.gz > ${sampleID}_raw_aln.sam

ml purge

### SAM-BAM
echo "Converting SAM-BAM"
ml samtools

samtools view -b ${sampleID}_raw_aln.sam > ${sampleID}_raw_aln.bam
samtools sort ${sampleID}_raw_aln.bam -o ${sampleID}_raw_aln_sorted.bam

ml purge


### Print time
end=`date +%s`
runtime=$((end-start))
out=$(echo "scale=2; (${runtime}/1)" | bc); echo "${out} seconds"
