#!/bin/bash
#SBATCH --job-name=pbcstat
#SBATCH --output=pbcstat_%j.out
#SBATCH --error=pbcstat_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

# Usage:
# sbatch pbcstat.slurm assembly.fasta ont_reads.fastq.gz outdir

ASSEMBLY="$1"
READS="$2"
OUTDIR="$3"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/purge_dups

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Map ONT reads to assembly
minimap2 \
    -x map-ont \
    -t "$THREADS" \
    "$ASSEMBLY" \
    "$READS" | gzip -c > ont.paf.gz

# Calculate base-level coverage statistics
pbcstat ont.paf.gz

echo "Done."