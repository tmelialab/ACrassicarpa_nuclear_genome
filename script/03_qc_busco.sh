#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=03_busco
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --time=24:00:00

set -euo pipefail

ENV_NAME="training_qc"
GENOME_FASTA="nd.asm.fasta"
OUT_NAME="hasil_busco"
LINEAGE="embryophyta_odb10"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mulai BUSCO (mode genome, lineage ${LINEAGE})..."
busco -i "${GENOME_FASTA}" \
      -o "${OUT_NAME}" \
      -m genome \
      -l "${LINEAGE}" \
      --cpu "${SLURM_CPUS_PER_TASK}" \
      --offline \
      --force

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Membuat visualisasi BUSCO plot..."
busco --plot "${OUT_NAME}"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] BUSCO selesai."
