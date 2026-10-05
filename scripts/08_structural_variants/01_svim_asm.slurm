#!/bin/bash
#SBATCH --job-name=svim_asm
#SBATCH --output=svim_asm_%j.out
#SBATCH --error=svim_asm_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

# Usage:
# sbatch run_svim_asm.slurm bam reference.fa outdir

BAM="$1"
REF="$2"
OUTDIR="$3"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate $HOME/.conda/envs/svim_asm

mkdir -p "$OUTDIR"

svim-asm haploid \
    --min_sv_size 50 \
    "$OUTDIR" \
    "$BAM" \
    "$REF"

echo "Done."
ls -lh "$OUTDIR"

# ============================================================
# Summarize and remove exact duplicate SV records
# ============================================================

VCF="${OUTDIR}/variants.vcf"
DEDUP_VCF="${OUTDIR}/variants.exact_dedup.vcf"
SUMMARY="${OUTDIR}/svim_summary.tsv"
STATS="${OUTDIR}/svim_summary.txt"

python3 - "$VCF" "$DEDUP_VCF" "$SUMMARY" "$STATS" <<'PY'
import sys
from collections import defaultdict

vcf, dedup_vcf, summary_file, stats_file = sys.argv[1:]


def parse_info(info_string):
    info = {}
    for item in info_string.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info


def get_length(pos, svtype, info):
    # Prefer SVLEN when available
    if "SVLEN" in info:
        try:
            return abs(int(info["SVLEN"].split(",")[0]))
        except ValueError:
            pass

    # Otherwise derive from END for interval-type SVs
    if "END" in info and svtype != "BND":
        try:
            return abs(int(info["END"]) - int(pos))
        except ValueError:
            pass

    return None


raw_count = defaultdict(int)
raw_bp = defaultdict(int)

dedup_count = defaultdict(int)
dedup_bp = defaultdict(int)

seen = set()

total_raw = 0
total_dedup = 0
duplicates = 0


with open(vcf) as infile, open(dedup_vcf, "w") as outfile:

    for line in infile:

        # Preserve VCF headers
        if line.startswith("#"):
            outfile.write(line)
            continue

        fields = line.rstrip("\n").split("\t")

        chrom = fields[0]
        pos   = fields[1]
        ref   = fields[3]
        alt   = fields[4]
        info_string = fields[7]

        info = parse_info(info_string)

        svtype = info.get("SVTYPE", "UNKNOWN")
        end    = info.get("END", "")
        svlen  = info.get("SVLEN", "")

        length = get_length(pos, svtype, info)

        # -----------------------
        # Raw statistics
        # -----------------------
        total_raw += 1
        raw_count[svtype] += 1

        if length is not None:
            raw_bp[svtype] += length

        # Exact duplicate definition:
        # same chromosome, position, REF, ALT,
        # SVTYPE, END and SVLEN
        key = (
            chrom,
            pos,
            ref,
            alt,
            svtype,
            end,
            svlen
        )

        if key in seen:
            duplicates += 1
            continue

        seen.add(key)

        outfile.write(line)

        total_dedup += 1
        dedup_count[svtype] += 1

        if length is not None:
            dedup_bp[svtype] += length


# ------------------------------------------------------------
# Write summary table
# ------------------------------------------------------------

all_types = sorted(set(raw_count) | set(dedup_count))

with open(summary_file, "w") as out:

    out.write(
        "SVTYPE\tRAW_COUNT\tRAW_TOTAL_BP\t"
        "DEDUP_COUNT\tDEDUP_TOTAL_BP\n"
    )

    for svtype in all_types:

        out.write(
            f"{svtype}\t"
            f"{raw_count[svtype]}\t"
            f"{raw_bp[svtype]}\t"
            f"{dedup_count[svtype]}\t"
            f"{dedup_bp[svtype]}\n"
        )


# ------------------------------------------------------------
# Overall statistics
# ------------------------------------------------------------

with open(stats_file, "w") as out:

    out.write(f"Raw SV records:             {total_raw}\n")
    out.write(f"Exact duplicate records:    {duplicates}\n")
    out.write(f"Nonredundant SV records:    {total_dedup}\n")

    if total_raw > 0:
        pct = 100 * duplicates / total_raw
        out.write(f"Exact duplicate percentage: {pct:.3f}%\n")


print("\n=== SVIM SUMMARY ===")
print(f"Raw SV records:          {total_raw}")
print(f"Exact duplicates:        {duplicates}")
print(f"Nonredundant SV records: {total_dedup}")

print("\nSummary:")
with open(summary_file) as f:
    print(f.read())
PY


echo "========================================="
echo "SVIM-asm analysis completed"
echo "========================================="
echo
cat "$STATS"
echo
cat "$SUMMARY"