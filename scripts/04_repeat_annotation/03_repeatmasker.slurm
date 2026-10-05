#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=06_repeatmasker
#SBATCH --output=logs/repeatmasker_%j.log
#SBATCH --error=logs/repeatmasker_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=short

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate repeats
cd "$(dirname "${BASH_SOURCE[0]:-$0}")"  # = script/ ; submit dari project root bila perlu

LIB="acacia_db-families.fa"
if [ ! -f "$LIB" ]; then
    LIB="$(ls RM_*/consensi.fa 2>/dev/null | head -1)"
fi
echo "--- RepeatMasker start $(date), lib=${LIB} ---"
RepeatMasker -pa 16 -lib "$LIB" -gff -xsmall acacia_final.fasta
echo "--- RepeatMasker done $(date) ---"
