#!/bin/bash
#SBATCH --job-name=del_repeat_enrich
#SBATCH --output=del_repeat_enrich_%j.out
#SBATCH --error=del_repeat_enrich_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

# Usage:
# sbatch del_repeat_enrichment.slurm \
#   del_repeat_enrichment.py \
#   supported_deletions.bed \
#   acra3rx_repeats.bed \
#   Acra_USDA_v1.100k.fa.fai \
#   output_dir \
#   [nperm]

SCRIPT="$1"
DELETIONS="$2"
REPEATS="$3"
FAI="$4"
OUTDIR="$5"
NPERM="${6:-1000}"

mkdir -p "$OUTDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/repeat

echo "Starting deletion-repeat enrichment"
echo "Permutations: $NPERM"
echo "Memory: ${SLURM_MEM_PER_NODE:-NA} MB"
echo

python3 "$SCRIPT" \
  --deletions "$DELETIONS" \
  --repeats "$REPEATS" \
  --fai "$FAI" \
  --outdir "$OUTDIR" \
  --nperm "$NPERM" \
  --seed 42

echo
echo "Finished."
echo

cat "$OUTDIR/deletion_repeat_enrichment.tsv"