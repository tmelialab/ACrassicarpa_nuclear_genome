#!/bin/bash
#SBATCH --job-name=get_seqs
#SBATCH --output=get_seqs_%j.out
#SBATCH --error=get_seqs_%j.err
#SBATCH --cpus-per-task=1


set -euo pipefail

# Usage:
# sbatch get_seqs.slurm dups.bed assembly.fasta outdir [prefix]

DUPS="$1"
ASSEMBLY="$2"
OUTDIR="$3"
PREFIX="${4:-acacia_purged}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/purge_dups

mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Conservative filtering:
# Remove HAPLOTIG, OVLP, and JUNK
# Keep REPEAT and HIGHCOV in the primary assembly
awk '$4=="HAPLOTIG" || $4=="OVLP" || $4=="JUNK"' \
    "$DUPS" > dups_filtered.bed

echo "Regions used for purging:"
cut -f4 dups_filtered.bed | sort | uniq -c

get_seqs \
    -e \
    -c \
    -p "$PREFIX" \
    dups_filtered.bed \
    "$ASSEMBLY"

echo "Done."
ls -lh