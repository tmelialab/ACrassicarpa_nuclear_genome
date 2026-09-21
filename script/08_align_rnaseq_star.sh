#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=08_star_align
#SBATCH --output=logs/star_align_%j.log
#SBATCH --error=logs/star_align_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=short

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate training_qc

GENOME="$HOME/rerunacacia/annotation/acacia_final.fasta.masked"
RNA_DIR="$HOME/rerunacacia/annotation/rnaseq"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

ulimit -n 65535
mkdir -p "${RNA_DIR}/star_index"

echo "--- 1. Generating STAR Genome Index $(date) ---"
STAR --runThreadN "$THREADS" \
     --runMode genomeGenerate \
     --genomeDir "${RNA_DIR}/star_index" \
     --genomeFastaFiles "$GENOME" \
     --genomeSAindexNbases 13

echo "--- 2. Aligning Paired-End Reads $(date) ---"
STAR --runThreadN "$THREADS" \
     --genomeDir "${RNA_DIR}/star_index" \
     --readFilesIn "${RNA_DIR}/SRR25080411_1.fastq" "${RNA_DIR}/SRR25080411_2.fastq" \
     --outFileNamePrefix "${RNA_DIR}/SRR25080411_" \
     --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 60000000000

STAR --runThreadN "$THREADS" \
     --genomeDir "${RNA_DIR}/star_index" \
     --readFilesIn "${RNA_DIR}/SRR25080412_1.fastq" "${RNA_DIR}/SRR25080412_2.fastq" \
     --outFileNamePrefix "${RNA_DIR}/SRR25080412_" \
     --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 60000000000

STAR --runThreadN "$THREADS" \
     --genomeDir "${RNA_DIR}/star_index" \
     --readFilesIn "${RNA_DIR}/SRR25816559_1.fastq" "${RNA_DIR}/SRR25816559_2.fastq" \
     --outFileNamePrefix "${RNA_DIR}/SRR25816559_" \
     --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 60000000000

echo "--- 3. Aligning Single-End Reads $(date) ---"
STAR --runThreadN "$THREADS" \
     --genomeDir "${RNA_DIR}/star_index" \
     --readFilesIn "${RNA_DIR}/SRR1168433.fastq" \
     --outFileNamePrefix "${RNA_DIR}/SRR1168433_" \
     --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM 60000000000

echo "--- Alignment Complete! $(date) ---"
