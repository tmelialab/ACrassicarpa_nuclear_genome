#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=11_trnascan
#SBATCH --output=logs/tRNAscan_%j.log
#SBATCH --error=logs/tRNAscan_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --partition=short

set -euo pipefail

ENV_NAME="training_qc"
GENOME_MASKED="$HOME/rerunacacia/annotation/acacia_final.fasta.masked"
OUTDIR="$HOME/rerunacacia/annotation/trnascan"

echo "=== Mulai tRNAscan-SE $(date) ==="
echo "Genome: $GENOME_MASKED ($(grep -c "^>" "$GENOME_MASKED") contigs)"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

command -v tRNAscan-SE >/dev/null || { echo "ERROR: tRNAscan-SE belum terinstall!"; exit 1; }
tRNAscan-SE -h 2>&1 | head -3 || true

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

echo "Menjalankan tRNAscan-SE (mode Eukaryotic) pada genome soft-masked..."
tRNAscan-SE \
    -E \
    -o trna.out \
    -m trna.stats \
    -f trna.fa \
    --thread "${SLURM_CPUS_PER_TASK}" \
    "${GENOME_MASKED}"

echo "=== tRNAscan-SE selesai $(date) ==="
echo "Output: $OUTDIR/trna.out ($(wc -l < trna.out) baris), trna.stats, trna.fa"
# ringkas jumlah tRNA
grep -c "tRNA" trna.out || true
cat trna.stats 2>&1 | head -20
