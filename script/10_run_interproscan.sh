#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=10_ipr_annot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --output=logs/ipr_annot_%j.log
#SBATCH --error=logs/ipr_annot_%j.err

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate training_qc

############################################
# USER SETTINGS
############################################

IPR=$HOME/acacia_project/tools/interproscan-5.75-106.0
INPUT=$HOME/rerunacacia/annotation/braker/braker.aa
OUTDIR=$HOME/rerunacacia/annotation/interproscan
PREFIX=braker
CPUS=${SLURM_CPUS_PER_TASK:-64}
TMPDIR_BASE="/mgpfs/scratch/$USER/tmp_ipr_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR_BASE"

############################################
# PREP
############################################

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "[$(date)] Starting InterProScan 64c"
echo "Input: $INPUT ($(grep -c "^>" "$INPUT" 2>/dev/null || echo ?) proteins)"
CLEANED="${INPUT%.aa}.clean.aa"
echo "[$(date)] Cleaning * from proteins (InterProScan cannot handle *)"
tr -d '*' < "$INPUT" > "$CLEANED"
INPUT="$CLEANED"
echo "Cleaned: $(grep -c "^>" "$INPUT") proteins, * removed"
echo "CPUs: $CPUS, TMPDIR: $TMPDIR_BASE"
java -version 2>&1 | head -1

############################################
# STEP 1 — RUN INTERPROSCAN
############################################

"$IPR"/interproscan.sh \
  -i "$INPUT" \
  -f tsv \
  -dp \
  -goterms \
  -iprlookup \
  -pa \
  -T "$TMPDIR_BASE" \
  -cpu "$CPUS" \
  -appl Pfam,SMART,ProSiteProfiles,ProSitePatterns,SUPERFAMILY \
  -b "$PREFIX"

echo "[$(date)] InterProScan finished, cleaning TMPDIR"
rm -rf "$TMPDIR_BASE"

############################################
# STEP 2 — SUMMARY STATISTICS
############################################

TSV="${PREFIX}.tsv"

echo "[$(date)] Generating summary statistics"

TOTAL_PROTEINS=$(grep -c "^>" "$INPUT")
PROTEINS_WITH_HITS=$(cut -f1 "$TSV" | sort | uniq | wc -l)
GO_TERMS=$(grep -o "GO:[0-9]\+" "$TSV" | sort | uniq | wc -l)
PROTEINS_WITH_GO=$(grep "GO:" "$TSV" | cut -f1 | sort | uniq | wc -l)
awk '$4=="Pfam"{print $6}' "$TSV" | sort | uniq -c | sort -nr  > pfam.txt
head -20 pfam.txt > top20_pfam.txt

############################################
# WRITE SUMMARY FILE
############################################

SUMMARY=interproscan_summary.txt

{
  echo "InterProScan summary"
  echo "===================="
  echo "Input protein file: $INPUT"
  echo "InterProScan version: 5.75-106.0"
  echo "CPUs: $CPUS (was 16, now 64 — 4x speedup)"
  echo ""
  echo "Total proteins: $TOTAL_PROTEINS"
  echo "Proteins with InterPro/Pfam hits: $PROTEINS_WITH_HITS"
  echo "Proteins with GO terms: $PROTEINS_WITH_GO"
  echo "Unique GO terms: $GO_TERMS"
  echo ""
  echo "Top 20 Pfam domains: see top20_pfam.txt"
} > "$SUMMARY"

echo "[$(date)] Summary written to $SUMMARY"
echo "[$(date)] Done."
