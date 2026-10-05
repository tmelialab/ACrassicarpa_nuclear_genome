#!/bin/bash -l
#SBATCH --job-name=nanoplot
#SBATCH --output=nanoplot_%j.out
#SBATCH --error=nanoplot_%j.err
#SBATCH --cpus-per-task=16

# AcraUNRI Oxford Nanopore read QC using NanoPlot
#
# Usage:
#   sbatch 01_nanoplot_ont.slurm <reads.fastq.gz> <output_directory>

set -euo pipefail

FASTQ="${1:?Usage: sbatch $0 <reads.fastq.gz> <output_directory>}"
OUTDIR="${2:?Usage: sbatch $0 <reads.fastq.gz> <output_directory>}"

echo "=========================================================="
echo "Start date : $(date)"
echo "Job name   : ${SLURM_JOB_NAME}"
echo "Job ID     : ${SLURM_JOB_ID}"
echo "Input      : ${FASTQ}"
echo "Output     : ${OUTDIR}"
echo "Threads    : ${SLURM_CPUS_PER_TASK}"
echo "=========================================================="

module load miniconda
conda activate nanoplot

mkdir -p "${OUTDIR}"

NanoPlot \
    --fastq "${FASTQ}" \
    --outdir "${OUTDIR}" \
    --threads "${SLURM_CPUS_PER_TASK}"

echo "=========================================================="
echo "End date : $(date)"
echo "=========================================================="
