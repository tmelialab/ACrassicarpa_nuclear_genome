#!/bin/bash
# NOTE (sanitasi Sep 2026): path absolut asal (home user asal,
# /mgpfs/scratch/$USER) diganti $HOME dan /mgpfs/scratch/$USER.
# Sesuaikan PROJECT_DIR / tool paths di bawah dengan akun HPC masing-masing.
#SBATCH --job-name=90_syri_plotsr
#SBATCH --output=syri_plotsr_%j.log
#SBATCH --error=syri_plotsr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=24:00:00

set -euo pipefail

echo "=================================================================="
echo "   SyRI + Visualisasi Perbandingan Genome"
echo "=================================================================="

# Load conda dan aktifkan environment syri
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate syri_env

cd $HOME/combine/compare

REF_FA="ref_acacia_crassicarpa.fa"
QRY_FA="my_assembly.fa"
ALN_BAM="alignment.sorted.bam"
THREADS="${SLURM_CPUS_PER_TASK:-32}"

for f in "${REF_FA}" "${QRY_FA}"; do
    test -s "${f}" || { echo "ERROR: File tidak ditemukan/kosong: ${f}"; exit 1; }
done

command -v minimap2 >/dev/null 2>&1 || { echo "ERROR: minimap2 tidak ditemukan di environment"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools tidak ditemukan di environment"; exit 1; }

echo "[3] Membuat alignment BAM dengan minimap2 --eqx..."
echo "Melakukan alignment dengan minimap2..."
minimap2 -ax asm5 --eqx -t "${THREADS}" "${REF_FA}" "${QRY_FA}" > alignment.sam
samtools view -bS alignment.sam | samtools sort -@ "${THREADS}" -o "${ALN_BAM}" -
samtools index "${ALN_BAM}"

echo "[3b] Menyusun pasangan kromosom 1:1 untuk SyRI..."
samtools view -F 2308 "${ALN_BAM}" | awk '
function refbp(cigar,   a,n,op,s) {
    s=0
    while (match(cigar, /([0-9]+)([MIDNSHP=X])/, a)) {
        n=a[1]+0
        op=a[2]
        if (op=="M" || op=="=" || op=="X" || op=="D" || op=="N") s+=n
        cigar=substr(cigar, RSTART+RLENGTH)
    }
    return s
}
$3 != "*" {
    key=$3 "\t" $1
    aln[key]+=refbp($6)
}
END {
    for (k in aln) print k "\t" aln[k]
}' | sort -k3,3nr > chr_pairs.raw.tsv

awk 'BEGIN{OFS="\t"}
{
    r=$1; q=$2; b=$3
    if (!(r in seenr) && !(q in seenq)) {
        print r, q, b
        seenr[r]=1
        seenq[q]=1
    }
}' chr_pairs.raw.tsv > chr_pairs.tsv

test -s chr_pairs.tsv || { echo "ERROR: Tidak ada pasangan kromosom 1:1 yang bisa dibuat untuk SyRI"; exit 1; }

REF_SYRI_FA="ref_for_syri.fa"
QRY_SYRI_FA="qry_for_syri.fa"
ALN_BAM_SYRI="alignment.syri.sorted.bam"
SYRI_PREFIX="syri_acacia_"

: > "${REF_SYRI_FA}"
: > "${QRY_SYRI_FA}"

while IFS=$'\t' read -r REF_ID QRY_ID _; do
    samtools faidx "${REF_FA}" "${REF_ID}" >> "${REF_SYRI_FA}"
    samtools faidx "${QRY_FA}" "${QRY_ID}" | awk -v rid="${REF_ID}" 'NR==1{$0=">"rid}1' >> "${QRY_SYRI_FA}"
done < chr_pairs.tsv

samtools faidx "${REF_SYRI_FA}"
samtools faidx "${QRY_SYRI_FA}"

echo "[3c] Realignment pada FASTA ter-normalisasi untuk SyRI..."
minimap2 -ax asm5 --eqx -t "${THREADS}" "${REF_SYRI_FA}" "${QRY_SYRI_FA}" | \
    samtools view -bS - | \
    samtools sort -@ "${THREADS}" -o "${ALN_BAM_SYRI}" -
samtools index "${ALN_BAM_SYRI}"

echo "[4] Menjalankan SyRI..."

export ALN_BAM_SYRI REF_SYRI_FA QRY_SYRI_FA SYRI_PREFIX
python - <<'PY'
import os
import sys
import pandas as pd

# Force pandas to disable copy-on-write before SyRI imports manipulate DataFrames.
pd.options.mode.copy_on_write = False

from syri.scripts.syri import main

sys.argv = [
    "syri",
    "-c", os.environ["ALN_BAM_SYRI"],
    "-F", "B",
    "--cigar",
    "--no-chrmatch",
    "--nc", "1",
    "-r", os.environ["REF_SYRI_FA"],
    "-q", os.environ["QRY_SYRI_FA"],
    "--prefix", os.environ["SYRI_PREFIX"],
]
main()
PY

echo "[4] SyRI selesai."

SYRI_OUT="${SYRI_PREFIX}syri.out"
test -s "${SYRI_OUT}"

echo "[5] Membuat visualisasi dengan plotsr..."

plotsr \
    --sr "${SYRI_OUT}" \
    --genomes "${QRY_SYRI_FA},${REF_SYRI_FA}" \
    -o acacia_comparison_plotsr.pdf \
    -W 14 \
    -H 10

echo "=================================================================="
echo "SELESAI!"
echo "Visualisasi: acacia_comparison_plotsr.pdf"
echo "Detail SV   : ${SYRI_OUT}"
echo "=================================================================="