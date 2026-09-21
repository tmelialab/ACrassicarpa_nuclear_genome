#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=07_sra_dump
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --time=24:00:00
#SBATCH --output=sra_dump_%j.out
#SBATCH --error=sra_dump_%j.err

set -euo pipefail

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Aktivasi Lingkungan Conda
: "${ENV_NAME:=training_qc}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
log_message "Conda environment '${ENV_NAME}' activated."

# Loop download dari sra_dump.txt
while read -r srr; do
    [ -z "$srr" ] && continue
    log_message "Downloading $srr..."
    fasterq-dump "$srr" --split-files --threads "${SLURM_CPUS_PER_TASK:-64}"
done < sra_dump.txt

log_message "Processing complete!"