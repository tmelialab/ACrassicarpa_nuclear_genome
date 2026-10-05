#!/bin/bash
#SBATCH --job-name=sv_gene_depletion
#SBATCH --output=sv_gene_depletion_%j.out
#SBATCH --error=sv_gene_depletion_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00

set -euo pipefail

SCRIPT="$1"
VCF="$2"
GTF="$3"
FAI="$4"
OUTDIR="$5"

mkdir -p "$OUTDIR"

python3 "$SCRIPT" \
  --vcf "$VCF" \
  --gtf "$GTF" \
  --fai "$FAI" \
  --outdir "$OUTDIR" \
  --nperm 1000 \
  --seed 42