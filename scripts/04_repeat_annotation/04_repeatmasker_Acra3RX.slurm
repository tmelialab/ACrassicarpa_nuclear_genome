#!/bin/bash
#SBATCH --job-name=repeatmasker_acra3rx
#SBATCH --output=repeatmasker_acra3rx_%j.out
#SBATCH --error=repeatmasker_acra3rx_%j.err
#SBATCH --cpus-per-task=8

# Repeat annotation of the Acra3RX reference using the
# species-specific A. crassicarpa repeat library generated
# from AcraUNRI with RepeatModeler.
#
# Usage:
#   sbatch 04_repeatmasker_Acra3RX.slurm \
#       <Acra3RX.fa> <repeat_library.fa> <output_directory>

set -euo pipefail

GENOME="${1:?Usage: sbatch $0 <Acra3RX.fa> <repeat_library.fa> <output_dir>}"
LIB="${2:?Usage: sbatch $0 <Acra3RX.fa> <repeat_library.fa> <output_dir>}"
OUTDIR="${3:?Usage: sbatch $0 <Acra3RX.fa> <repeat_library.fa> <output_dir>}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate repeat
fi

command -v RepeatMasker >/dev/null 2>&1 || {
    echo "ERROR: RepeatMasker is not available in the current environment." >&2
    exit 1
}

mkdir -p "${OUTDIR}"

echo "=========================================================="
echo "Genome         : ${GENOME}"
echo "Repeat library : ${LIB}"
echo "Output         : ${OUTDIR}"
echo "Threads        : ${THREADS}"
echo "Start          : $(date)"
echo "=========================================================="

RepeatMasker \
    -pa "${THREADS}" \
    -gff \
    -lib "${LIB}" \
    -dir "${OUTDIR}" \
    -xsmall \
    "${GENOME}"

echo "=========================================================="
echo "Finished: $(date)"
echo "=========================================================="
