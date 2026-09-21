#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=13_diamond_uniref
#SBATCH --output=logs/diamond_uniref_%j.log
#SBATCH --error=logs/diamond_uniref_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=20:00:00
#SBATCH --partition=short

set -euo pipefail

ENV_NAME="training_qc"
QUERY="$HOME/rerunacacia/annotation/braker/braker.clean.aa"
if [ ! -f "$QUERY" ]; then
  echo "Cleaning * from braker.aa -> $QUERY"
  tr -d '*' < $HOME/rerunacacia/annotation/braker/braker.aa > "$QUERY"
fi
DB="$HOME/acacia_project/refdata/uniref90/uniref90db.dmnd"
if [ ! -f "$DB" ]; then
  DB="$HOME/acacia_project/refdata/uniref90.fasta"
fi
OUTDIR="$HOME/rerunacacia/annotation/diamond"
OUT_TSV="${OUTDIR}/diamond_uniref90.tsv"

echo "=== Mulai Diamond BLASTp vs UniRef90 $(date) ==="
echo "Query: $QUERY ($(grep -c "^>" "$QUERY") proteins)"
echo "DB: $DB ($(ls -lh "$DB" | awk "{print \$5}"))"
echo "Out: $OUT_TSV"
mkdir -p "$OUTDIR"
ls -lh "$DB" | head -1

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"
diamond --version 2>&1 | head -1

diamond blastp \
    --query "${QUERY}" \
    --db "${DB}" \
    --out "${OUT_TSV}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive \
    --max-target-seqs 1 \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --evalue 1e-5

echo "=== Diamond BLASTp UniRef90 selesai $(date) ==="
echo "Jumlah protein teranotasi:"
cut -f1 "${OUT_TSV}" | sort | uniq | wc -l
echo "Hasil: ${OUT_TSV} ($(wc -l < "${OUT_TSV}") hits)"
ls -lh "${OUT_TSV}"
