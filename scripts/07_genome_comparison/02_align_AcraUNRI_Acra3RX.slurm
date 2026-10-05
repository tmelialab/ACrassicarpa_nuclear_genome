#!/bin/bash
#SBATCH --job-name=acra3rx_align
#SBATCH --output=acra3rx_align_%j.out
#SBATCH --error=acra3rx_align_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

# Usage:
# sbatch align_acra3rx.slurm reference.fa query.fa outdir [prefix]

REF="$1"
QUERY="$2"
OUTDIR="$3"
PREFIX="${4:-acra_vs_acra3rx}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

module load bioinformatics/minimap2/2.28

mkdir -p "$OUTDIR"

BAM="${OUTDIR}/${PREFIX}.sorted.bam"

minimap2 \
    -a \
    -x asm5 \
    --cs \
    -r2k \
    -t "$THREADS" \
    "$REF" \
    "$QUERY" | \
samtools sort \
    -@ "$THREADS" \
    -m 3G \
    -o "$BAM"

samtools index "$BAM"

samtools flagstat "$BAM" \
    > "${OUTDIR}/${PREFIX}.flagstat.txt"

samtools coverage "$BAM" \
    > "${OUTDIR}/${PREFIX}.coverage.tsv"

echo "Done."