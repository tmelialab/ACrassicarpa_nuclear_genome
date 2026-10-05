#!/bin/bash
#SBATCH --job-name=busco_proteins
#SBATCH --output=busco_proteins_%j.out
#SBATCH --error=busco_proteins_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# BUSCO completeness assessment of the representative
# AcraUNRI predicted protein set.
#
# Usage:
#   sbatch 04_busco_proteins.slurm \
#       <proteins.fa> <embryophyta_odb10_directory> <output_directory>

set -euo pipefail

PROTEINS="${1:?Usage: sbatch $0 <proteins.fa> <lineage_dir> <output_dir>}"
LINEAGE="${2:?Usage: sbatch $0 <proteins.fa> <lineage_dir> <output_dir>}"
OUTDIR="${3:?Usage: sbatch $0 <proteins.fa> <lineage_dir> <output_dir>}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Protein set       : ${PROTEINS}"
echo "BUSCO lineage     : ${LINEAGE}"
echo "Output directory  : ${OUTDIR}"
echo "Threads           : ${THREADS}"
echo "=========================================================="

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate busco
fi

command -v busco >/dev/null 2>&1 || {
    echo "ERROR: BUSCO is not available in the current environment." >&2
    exit 1
}

mkdir -p "${OUTDIR}"

busco \
    -i "${PROTEINS}" \
    -l "${LINEAGE}" \
    -m protein \
    -o busco_proteins_acraunri \
    -c "${THREADS}" \
    --out_path "${OUTDIR}" \
    --offline \
    -f

echo "=========================================================="
echo "Finished on: $(date)"
echo "=========================================================="
