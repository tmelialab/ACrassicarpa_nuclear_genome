#!/bin/bash
#SBATCH --job-name=meryl_kmer
#SBATCH --output=meryl_kmer_%j.out
#SBATCH --error=meryl_kmer_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

# Usage:
# sbatch 01_meryl_kmer.sh R1.fastq.gz R2.fastq.gz output_directory [k]
#
# Example:
# sbatch 01_meryl_kmer.sh \
#   /mgpfs/home/tishalia/acra/data/ACGenome_ILM/Raw_97_R1.fq.gz \
#   /mgpfs/home/tishalia/acra/data/ACGenome_ILM/Raw_97_R2.fq.gz \
#   /mgpfs/home/tishalia/acra/kmer_analysis \
#   21

R1="$1"
R2="$2"
OUTDIR="$3"
K="${4:-21}"

THREADS="${SLURM_CPUS_PER_TASK:-32}"

mkdir -p "${OUTDIR}"

# Activate environment containing meryl
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kmer

echo "[$(date)] Counting ${K}-mers..."

meryl \
    k="${K}" \
    threads="${THREADS}" \
    count \
    output "${OUTDIR}/reads_k${K}.meryl" \
    "${R1}" "${R2}"

echo "[$(date)] Generating histogram..."

meryl histogram \
    "${OUTDIR}/reads_k${K}.meryl" \
    > "${OUTDIR}/reads_k${K}.hist"

echo "[$(date)] Done."
echo "Histogram: ${OUTDIR}/reads_k${K}.hist"