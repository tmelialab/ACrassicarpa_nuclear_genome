#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=12_diamond_sp
#SBATCH --output=logs/diamond_sp_%j.log
#SBATCH --error=logs/diamond_sp_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=short

set -euo pipefail

ENV_NAME="training_qc"
QUERY="$HOME/rerunacacia/annotation/braker/braker.clean.aa"
# fallback: if clean not exists, create from braker.aa
if [ ! -f "$QUERY" ]; then
  echo "Cleaning * from braker.aa -> $QUERY"
  tr -d '*' < $HOME/rerunacacia/annotation/braker/braker.aa > "$QUERY"
fi
DB="$HOME/acacia_project/refdata/swissprot/swissprot.dmnd"
# alternative DB path if above not found: check uniprot_sprot
if [ ! -f "$DB" ]; then
  DB="$HOME/acacia_project/refdata/swissprot.dmnd"
fi
OUTDIR="$HOME/rerunacacia/annotation/diamond"
OUT_TSV="${OUTDIR}/diamond_swissprot.tsv"

echo "=== Mulai Diamond BLASTp vs Swiss-Prot $(date) ==="
echo "Query: $QUERY ($(grep -c "^>" "$QUERY") proteins)"
echo "DB: $DB ($(ls -lh "$DB" | awk "{print \$5}"))"
echo "Out: $OUT_TSV"
mkdir -p "$OUTDIR"

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

echo "=== Diamond BLASTp Swiss-Prot selesai $(date) ==="
echo "Jumlah protein teranotasi:"
cut -f1 "${OUT_TSV}" | sort | uniq | wc -l
echo "Hasil: ${OUT_TSV} ($(wc -l < "${OUT_TSV}") hits)"
ls -lh "${OUT_TSV}"
