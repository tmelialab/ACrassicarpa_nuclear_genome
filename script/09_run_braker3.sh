#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=09_braker3
#SBATCH --output=logs/braker3_%j.log
#SBATCH --error=logs/braker3_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=short

set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate training_qc
export GENEMARK_PATH="$HOME/acacia_project/tools/etpbraker/bin"
export PROTHINT_PATH="$HOME/acacia_project/tools/prothint/ProtHint-2.6.0/bin"

GENOME="$HOME/rerunacacia/annotation/acacia_final.fasta.masked"
RNA_DIR="$HOME/rerunacacia/annotation/rnaseq"
PROT="$HOME/acacia_project/annotation/repeatmodeler/Viridiplantae_dwipa.fa"
WORKDIR="$HOME/rerunacacia/annotation/braker"
BAMS="${RNA_DIR}/SRR1168433_Aligned.sortedByCoord.out.bam,${RNA_DIR}/SRR25080411_Aligned.sortedByCoord.out.bam,${RNA_DIR}/SRR25080412_Aligned.sortedByCoord.out.bam,${RNA_DIR}/SRR25816559_Aligned.sortedByCoord.out.bam"

echo "--- BRAKER3 start $(date) ---"
echo "GENEMARK_PATH=$GENEMARK_PATH"
echo "PROTHINT_PATH=$PROTHINT_PATH"
ls -lh "$GENEMARK_PATH/gmetp.pl" | head -1
ls -lh "$PROTHINT_PATH/prothint.py" | head -1
braker.pl \
  --genome="$GENOME" \
  --bam="$BAMS" \
  --prot_seq="$PROT" \
  --softmasking \
  --threads=32 \
  --species=acacia_crassicarpa \
  --workingdir="$WORKDIR" \
  --gff3
echo "--- BRAKER3 done $(date) ---"
