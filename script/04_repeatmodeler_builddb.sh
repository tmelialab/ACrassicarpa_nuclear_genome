#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=04_builddb
#SBATCH --output=logs/builddb_%j.log
#SBATCH --error=logs/builddb_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=medium-small

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate repeats
cd "$(dirname "${BASH_SOURCE[0]:-$0}")"  # = script/ ; submit dari project root bila perlu

echo "--- BuildDatabase start $(date) ---"
BuildDatabase -name acacia_db  acacia_final.fasta
echo "--- BuildDatabase done $(date) ---"
