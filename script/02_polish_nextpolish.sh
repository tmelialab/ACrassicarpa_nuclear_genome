#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=02_nextpolish
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=240G
#SBATCH --time=23:00:00
#SBATCH --output=polish_%j.out
#SBATCH --error=polish_%j.err

set -euo pipefail

# NextPolish butuh libbz2.so.1.0 dari conda env di LD_LIBRARY_PATH
export PATH="$HOME/.conda/envs/training_qc/bin:$PATH"
export LD_LIBRARY_PATH="$HOME/.conda/envs/training_qc/lib:${LD_LIBRARY_PATH:-}"
export NEXT_POLISH="$HOME/acacia_project/tools/NextPolish"

CONFIG_FILE="02_polish.cfg"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Mulai NextPolish..."
${NEXT_POLISH}/nextPolish "${CONFIG_FILE}"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] NextPolish selesai."
