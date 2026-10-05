#!/bin/bash
#SBATCH --job-name=genomescope2
#SBATCH --output=genomescope2_%j.out
#SBATCH --error=genomescope2_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00

set -euo pipefail

# Usage:
# sbatch 03_genomescope2.sh reads_k21.hist output_directory [k] [ploidy]

HIST="$1"
OUTDIR="$2"
K="${3:-21}"
PLOIDY="${4:-2}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate kmer

mkdir -p "$OUTDIR"

genomescope2 \
    -i "$HIST" \
    -o "$OUTDIR" \
    -k "$K" \
    -p "$PLOIDY"

echo "Done."