#!/bin/bash
#SBATCH --job-name=purgedups_self
#SBATCH --output=purgedups_self_%j.out
#SBATCH --error=purgedups_self_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00

set -euo pipefail

# Usage:
# sbatch self_align.slurm assembly.fasta outdir

ASSEMBLY="$1"
OUTDIR="$2"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/purge_dups

mkdir -p "$OUTDIR"

SPLIT="${OUTDIR}/assembly.split.fa"
PAF="${OUTDIR}/assembly.self.paf.gz"

split_fa "$ASSEMBLY" > "$SPLIT"

minimap2 \
    -xasm5 \
    -DP \
    -I 16G \
    -t "$THREADS" \
    "$SPLIT" \
    "$SPLIT" | gzip -c > "$PAF"

echo "Done."