#!/bin/bash
#SBATCH --job-name=purge_dups
#SBATCH --output=purge_dups_%j.out
#SBATCH --error=purge_dups_%j.err
#SBATCH --cpus-per-task=1


set -euo pipefail

# Usage:
# sbatch 08_run_purge_dups.sh \
#   self_alignment.paf.gz \
#   PB.base.cov \
#   cutoffs \
#   outdir

PAF="$1"
BASECOV="$2"
CUTOFFS="$3"
OUTDIR="$4"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/purge_dups

mkdir -p "$OUTDIR"

purge_dups \
    -2 \
    -T "$CUTOFFS" \
    -c "$BASECOV" \
    "$PAF" \
    > "$OUTDIR/dups.bed" \
    2> "$OUTDIR/purge_dups.log"

echo "Done."

