#!/bin/bash
#SBATCH --job-name=sniffles_force
#SBATCH --output=sniffles_force_%j.out
#SBATCH --error=sniffles_force_%j.err
#SBATCH --cpus-per-task=8


set -euo pipefail

# Usage:
# sbatch sniffles_force.slurm BAM SVIM_VCF REF OUTDIR [PREFIX]

BAM="$1"
SVIM_VCF="$2"
REF="$3"
OUTDIR="$4"
PREFIX="${5:-acra_svim_force}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$OUTDIR"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/sniffles2

OUTVCF="${OUTDIR}/${PREFIX}.vcf.gz"

echo "Force-genotyping SVIM-asm SVs with ONT reads..."

sniffles \
    --input "$BAM" \
    --genotype-vcf "$SVIM_VCF" \
    --vcf "$OUTVCF" \
    --reference "$REF" \
    --threads "$THREADS" \
    --allow-overwrite

echo
echo "Force genotyping completed."
echo "Output: $OUTVCF"

echo
echo "Total variants:"
zgrep -vc '^#' "$OUTVCF"

echo
echo "Genotype counts:"
zgrep -v '^#' "$OUTVCF" | \
awk -F'\t' '
{
    split($10,a,":")
    gt[a[1]]++
}
END {
    for(g in gt)
        print g, gt[g]
}' | sort



# ============================================================
# Create ONT-supported high-confidence SV set
# ============================================================

FORCE_VCF="${OUTDIR}/${PREFIX}.vcf.gz"
HC_VCF="${OUTDIR}/${PREFIX}.supported.vcf"
HC_SUMMARY="${OUTDIR}/${PREFIX}.supported_summary.tsv"

python3 - "$FORCE_VCF" "$HC_VCF" "$HC_SUMMARY" <<'PY'
import sys
import gzip
from collections import defaultdict

input_vcf, output_vcf, summary_file = sys.argv[1:]

def parse_info(info_string):
    info = {}
    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info

counts = defaultdict(int)
gt_counts = defaultdict(int)

total = 0
supported = 0

with gzip.open(input_vcf, "rt") as inp, open(output_vcf, "w") as out:

    for line in inp:

        if line.startswith("#"):
            out.write(line)
            continue

        fields = line.rstrip("\n").split("\t")

        total += 1

        fmt = fields[8].split(":")
        sample = fields[9].split(":")

        fmt_dict = dict(zip(fmt, sample))
        gt = fmt_dict.get("GT", "./.")

        # Normalize phased genotypes
        gt_norm = gt.replace("|", "/")

        # Keep variants supported by the ONT reads
        if gt_norm not in {"0/1", "1/0", "1/1"}:
            continue

        supported += 1

        info = parse_info(fields[7])
        svtype = info.get("SVTYPE", "UNKNOWN")

        counts[svtype] += 1

        if gt_norm in {"0/1", "1/0"}:
            gt_counts["HET"] += 1
        elif gt_norm == "1/1":
            gt_counts["HOM_ALT"] += 1

        out.write(line)


with open(summary_file, "w") as out:

    out.write("SVTYPE\tSUPPORTED_COUNT\n")

    for svtype in sorted(counts):
        out.write(f"{svtype}\t{counts[svtype]}\n")

    out.write("\n")
    out.write(f"Total force-genotyped\t{total}\n")
    out.write(f"ONT-supported\t{supported}\n")
    out.write(f"Heterozygous\t{gt_counts['HET']}\n")
    out.write(f"Homozygous_alt\t{gt_counts['HOM_ALT']}\n")

    if total:
        out.write(
            f"Supported_percent\t{100*supported/total:.2f}\n"
        )

print("Total force-genotyped:", total)
print("ONT-supported:", supported)
print("Heterozygous:", gt_counts["HET"])
print("Homozygous-alt:", gt_counts["HOM_ALT"])

print("\nSupported SV types:")
for svtype in sorted(counts):
    print(svtype, counts[svtype])
PY

echo
echo "======================================"
echo "High-confidence SV summary"
echo "======================================"

cat "$HC_SUMMARY"