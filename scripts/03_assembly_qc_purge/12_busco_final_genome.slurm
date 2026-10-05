#!/bin/bash
#SBATCH --job-name=busco_acra_final
#SBATCH --output=busco_acra_final_%j.out
#SBATCH --error=busco_acra_final_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# BUSCO assessment of the final AcraUNRI genome assembly
#
# Usage:
#   sbatch 12_busco_final_genome.slurm \
#       <genome.fa> <embryophyta_odb10_directory> <output_directory>
#
# Example:
#   sbatch 12_busco_final_genome.slurm \
#       AcraUNRI_v1.fa \
#       busco_downloads/lineages/embryophyta_odb10 \
#       busco_final

set -euo pipefail

GENOME="${1:?Usage: sbatch $0 <genome.fa> <lineage_dir> <output_dir>}"
LINEAGE="${2:?Usage: sbatch $0 <genome.fa> <lineage_dir> <output_dir>}"
OUTDIR="${3:?Usage: sbatch $0 <genome.fa> <lineage_dir> <output_dir>}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "=========================================================="
echo "Starting on       : $(date)"
echo "Running on node   : $(hostname)"
echo "Genome            : ${GENOME}"
echo "BUSCO lineage     : ${LINEAGE}"
echo "Output directory  : ${OUTDIR}"
echo "Threads           : ${THREADS}"
echo "=========================================================="

# Activate the BUSCO environment if Conda is available.
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
    -i "${GENOME}" \
    -o busco_embryophyta_final \
    --out_path "${OUTDIR}" \
    -m genome \
    -l "${LINEAGE}" \
    -c "${THREADS}" \
    --offline \
    -f

echo "=========================================================="
echo "Finished on: $(date)"
echo "=========================================================="
