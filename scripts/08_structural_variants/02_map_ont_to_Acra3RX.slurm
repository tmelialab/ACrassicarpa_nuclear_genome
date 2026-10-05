#!/bin/bash
#SBATCH --job-name=ont_map
#SBATCH --output=ont_map_%j.out
#SBATCH --error=ont_map_%j.err
#SBATCH --cpus-per-task=16

set -euo pipefail

# Usage:
# sbatch map_ont.slurm genome.fasta ont_reads.fastq.gz outdir

GENOME="$1"
READS="$2"
OUTDIR="$3"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

module load bioinformatics/minimap2/2.28
mkdir -p "$OUTDIR"

PREFIX="${OUTDIR}/ont_vs_assembly"

# Map ONT reads
minimap2 -ax map-ont -t "$THREADS" "$GENOME" "$READS" | \
    samtools sort -@ "$THREADS" -o "${PREFIX}.bam"

samtools index "${PREFIX}.bam"

# Mapping statistics
samtools flagstat "${PREFIX}.bam" \
    > "${PREFIX}.flagstat.txt"

# Primary alignments only
samtools view -b -F 2308 "${PREFIX}.bam" \
    > "${PREFIX}.primary.bam"

samtools index "${PREFIX}.primary.bam"

# Per-contig coverage
samtools coverage "${PREFIX}.primary.bam" \
    > "${PREFIX}.coverage.tsv"

echo "Done."