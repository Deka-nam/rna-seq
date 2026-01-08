#!/bin/bash

#SBATCH --job-name=trim_KK                    # Job name
#SBATCH --output=trim_KK_%j.log               # Output log file (%j will be replaced by the job ID)
#SBATCH --error=trim_KK_%j.err                # Error log file
#SBATCH --time=01:00:00                        # Time limit hrs:min:sec (adjust as necessary)
#SBATCH --ntasks=1                             # Number of tasks (processes)
#SBATCH --cpus-per-task=16                      # Number of CPU cores per task
#SBATCH --mem=32G                               # Total memory limit
#SBATCH --partition=general-compute            # Partition name
#SBATCH --qos=general-compute                   # QoS name (change to a valid QoS)
#SBATCH --cluster=ub-hpc

# load trimmomatic
module load trimmomatic/0.39-Java-11.0.16

FASTQDIR="/path/to/fastqfiles"

# Iterate over R1 files
for r1 in ${FASTQDIR}/*_1.fastq; do
    # Define R2 by substitution
    r2="${r1/_1.fastq/_2.fastq}"

    # Get base sample name
    sample=$(basename "$r1" _1.fastq)

    # Run Trimmomatic
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
    -phred33 -threads 8 \
    "$r1" "$r2" \
    "${FASTQDIR}/${sample}_1.trim.fastq" \
    "${FASTQDIR}/${sample}_1un.trim.fastq" \
    "${FASTQDIR}/${sample}_2.trim.fastq" \
    "${FASTQDIR}/${sample}_2un.trim.fastq" \
    ILLUMINACLIP:combined_adapters.fasta:2:30:10 \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:70

    echo "Trimming complete for ${sample}"
done
