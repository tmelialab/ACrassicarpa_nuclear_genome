#!/bin/bash -l
#SBATCH --job-name=jellyfish_kmers
#SBATCH --output=jellyfish_kmers_%j.out
#SBATCH --error=jellyfish_kmers_%j.err
#SBATCH --cpus-per-task=16

# Generate a 21-mer histogram from AcraUNRI paired-end short reads
# for GenomeScope analysis.
#
# Usage:
#   sbatch 03_jellyfish_kmers.slurm \
#       <R1.fastq.gz> <R2.fastq.gz> <output_prefix>

set -euo pipefail

READ1="${1:?Usage: sbatch $0 <R1.fastq.gz> <R2.fastq.gz> <output_prefix>}"
READ2="${2:?Usage: sbatch $0 <R1.fastq.gz> <R2.fastq.gz> <output_prefix>}"
PREFIX="${3:?Usage: sbatch $0 <R1.fastq.gz> <R2.fastq.gz> <output_prefix>}"

KMER=21
HASH_SIZE=16000000000

OUT_JF="${PREFIX}.jf"
OUT_HIST="${PREFIX}.histo"

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : ${SLURM_JOB_NAME}"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Read 1     : ${READ1}"
echo "Read 2     : ${READ2}"
echo "K-mer size : ${KMER}"
echo "Threads    : ${SLURM_CPUS_PER_TASK}"
echo "=========================================================="

module load jellyfish

zcat "${READ1}" "${READ2}" | \
    jellyfish count \
        -C \
        -m "${KMER}" \
        -s "${HASH_SIZE}" \
        -t "${SLURM_CPUS_PER_TASK}" \
        -o "${OUT_JF}" \
        /dev/fd/0

jellyfish histo \
    -t "${SLURM_CPUS_PER_TASK}" \
    "${OUT_JF}" \
    > "${OUT_HIST}"

echo "=========================================================="
echo "Jellyfish database : ${OUT_JF}"
echo "K-mer histogram    : ${OUT_HIST}"
echo "End date            : $(date)"
echo "=========================================================="
