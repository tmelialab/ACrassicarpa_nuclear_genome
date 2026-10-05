#!/bin/bash
#SBATCH --job-name=calcuts
#SBATCH --output=calcuts_%j.out
#SBATCH --error=calcuts_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

set -euo pipefail

# Usage:
# sbatch calcuts.slurm PB.stat outdir

PBSTAT="$1"
OUTDIR="$2"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/purge_dups

mkdir -p "$OUTDIR"

calcuts "$PBSTAT" \
    > "$OUTDIR/cutoffs" \
    2> "$OUTDIR/calcuts.log"