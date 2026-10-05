#!/bin/bash
#SBATCH --job-name=merqury
#SBATCH --output=merqury_%j.out
#SBATCH --error=merqury_%j.err
#SBATCH --cpus-per-task=4

set -euo pipefail

# Usage:
# sbatch merqury.sh reads_k21.meryl assembly.fasta outdir prefix

MERYL_DB="$1"
ASSEMBLY="$2"
OUTDIR="$3"
PREFIX="$4"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/kmer

mkdir -p "$OUTDIR"
cd "$OUTDIR"

merqury.sh \
    "$MERYL_DB" \
    "$ASSEMBLY" \
    "$PREFIX"

echo "Done."