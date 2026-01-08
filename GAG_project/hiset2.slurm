#!/bin/bash

#SBATCH --job-name=hiset2        # Job name
#SBATCH --output=hiset2_%j.log   # Output log file (%j will be replaced by the job ID)
#SBATCH --error=hiset2_%j.err    # Error log file
#SBATCH --time=50:00:00                # Time limit hrs:min:sec (adjust as necessary)
#SBATCH --ntasks=1                     # Number of tasks (processes)
#SBATCH --cpus-per-task=32              # Number of CPU cores per task
#SBATCH --mem=32G                       # Total memory limit
#SBATCH --partition=scavenger    # Partition name
#SBATCH --qos=scavenger                 # QoS name (change to a valid QoS)
#SBATCH --cluster=ub-hpc


set -euo pipefail

# Load modules 
module load gcc
module load openmpi
module load hisat2
module load samtools

# Paths
INDEX="hisat_index/mm10/mm10_index"
READ_DIR="trimmed_pairedreads"
OUT_DIR="aligned_mm10"

mkdir -p "${OUT_DIR}"

shopt -s nullglob

READ_DIR=trimmed_pairedreads

for R1 in "${READ_DIR}"/*_1.trim.fastq; do
    SAMPLE=$(basename "$R1" _1.trim.fastq)
    R2="${READ_DIR}/${SAMPLE}_2.trim.fastq"

    echo "R1: $R1"
    echo "R2: $R2"

    [[ -f "$R2" ]] || { echo "Missing R2 for $SAMPLE"; continue; }

    echo "Processing $SAMPLE"
    # hisat2 + samtools here
done

