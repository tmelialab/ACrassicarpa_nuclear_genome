#!/bin/bash
#SBATCH --job-name=sv_gene_impact
#SBATCH --output=sv_gene_impact_%j.out
#SBATCH --error=sv_gene_impact_%j.err
#SBATCH --cpus-per-task=4

set -euo pipefail

# Usage:
# sbatch sv_gene_impact.slurm \
#   sv_gene_impact.py \
#   supported.vcf \
#   annotation.gtf \
#   output_dir \
#   [distance]

SCRIPT="$1"
VCF="$2"
GTF="$3"
OUTDIR="$4"
DISTANCE="${5:-5000}"

mkdir -p "$OUTDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/repeat

echo "========================================"
echo "SV-gene impact analysis"
echo "========================================"
echo "VCF:      $VCF"
echo "GTF:      $GTF"
echo "Distance: $DISTANCE bp"
echo "Output:   $OUTDIR"
echo

python3 "$SCRIPT" \
  --vcf "$VCF" \
  --gtf "$GTF" \
  --outdir "$OUTDIR" \
  --distance "$DISTANCE"

echo
echo "========================================"
echo "Analysis completed"
echo "========================================"

echo
cat "$OUTDIR/sv_gene_overview.txt"

echo
cat "$OUTDIR/sv_gene_impact_summary.tsv"