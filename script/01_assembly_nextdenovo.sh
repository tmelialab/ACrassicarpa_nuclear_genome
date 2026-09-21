#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=01_nextdenovo
#SBATCH --output=nextdenovo_%j.out
#SBATCH --error=nextdenovo_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G

set -euo pipefail

#--- Konfigurasi ---
CONFIG_FILE="01_nextdenovo.cfg"
ENV_NAME="training_qc"
INPUT_FASTQ="$HOME/merged.fastq.gz"
NP_DIR="nanoplot_out"

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

log_message "================ JOB START ================"
log_message "Job ID: ${SLURM_JOB_ID}"
log_message "Host: $(hostname)"
log_message "-------------------------------------------"

# Aktivasi Conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"
log_message "Conda environment '${ENV_NAME}' aktif."

# STEP 0: QC Pra-assembly dengan NanoPlot
if [ ! -s "${NP_DIR}/NanoStats.txt" ]; then
    log_message "Menjalankan QC pra-assembly (NanoPlot)..."
    mkdir -p "${NP_DIR}"
    NanoPlot -t "${SLURM_CPUS_PER_TASK}" \
             --fastq "${INPUT_FASTQ}" \
             -o "${NP_DIR}" \
             --maxlength 50000 \
             --plots hex dot
    log_message "NanoPlot selesai: ${NP_DIR}"
else
    log_message "Hasil NanoPlot ditemukan di ${NP_DIR}, lanjut ke NextDenovo."
fi

# STEP 1: NextDenovo Assembly
command -v nextDenovo >/dev/null || { echo "ERROR: nextDenovo tidak ditemukan di env '$ENV_NAME'"; exit 127; }

log_message "Memulai NextDenovo assembly..."
nextDenovo "${CONFIG_FILE}"

log_message "Assembly selesai. File output: hasil/03.ctg_graph/nd.asm.fasta"
log_message "================= JOB END ================="
