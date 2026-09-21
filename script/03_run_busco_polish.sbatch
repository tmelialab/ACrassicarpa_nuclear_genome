#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=busco_polish
#SBATCH --output=logs/busco_polish_%j.out
#SBATCH --error=logs/busco_polish_%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G

BASE_DIR="$HOME/rerunacacia"
ASSEMBLY="${BASE_DIR}/polish/round_5/01_rundir/genome.nextpolish.fasta"
LINEAGE="embryophyta_odb10"
OUT="hasil_polish_round5"

echo "[$(date +%F\ %T)] JOB START"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate training_qc
echo "[$(date +%F\ %T)] conda training_qc activated"

if [ ! -f "${ASSEMBLY}" ]; then
    echo "ERROR: assembly not found: ${ASSEMBLY}"
    exit 1
fi

cd "${BASE_DIR}/busco"
busco -i "${ASSEMBLY}" \
      -o "${OUT}" \
      -m genome \
      -l "${LINEAGE}" \
      --cpu "${SLURM_CPUS_PER_TASK}" \
      --offline \
      --force \
      --download_path "${BASE_DIR}/busco/busco_downloads"
EXIT=$?

echo "[$(date +%F\ %T)] BUSCO exit code: ${EXIT}"
busco --plot "${OUT}" || true

if [ ${EXIT} -eq 0 ]; then
    echo "SUCCESS: hasil di ${BASE_DIR}/busco/${OUT}"
else
    echo "ERROR: BUSCO failed exit ${EXIT}"
fi
echo "[$(date +%F\ %T)] JOB END"
