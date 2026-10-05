#!/bin/bash
#SBATCH --job-name=sv_repeat
#SBATCH --output=sv_repeat_%j.out
#SBATCH --error=sv_repeat_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

# Usage:
# sbatch sv_repeat_analysis.slurm supported.vcf Acra3RX.fa.out consensi.fa.classified output_dir

VCF="$1"
RM_REF="$2"
LIB="$3"
OUTDIR="$4"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTDIR" "$OUTDIR/insertion_repeatmasker"

command -v python3 >/dev/null || { echo "python3 not found"; exit 1; }
command -v bedtools >/dev/null || { echo "bedtools not found"; exit 1; }
command -v RepeatMasker >/dev/null || { echo "RepeatMasker not found"; exit 1; }

DEL_BED="$OUTDIR/supported_deletions.bed"
INS_FA="$OUTDIR/supported_insertions.fa"
INS_INFO="$OUTDIR/supported_insertions.tsv"
REF_RM_BED="$OUTDIR/acra3rx_repeats.bed"
DEL_INT="$OUTDIR/deletion_repeat_intersections.tsv"

# ============================================================
# 1. Extract supported DEL intervals and INS sequences
# ============================================================

python3 - "$VCF" "$DEL_BED" "$INS_FA" "$INS_INFO" <<'PY'
import sys, gzip

vcf, delbed, insfa, insinfo = sys.argv[1:]

def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def info_dict(s):
    d = {}
    for x in s.split(";"):
        if "=" in x:
            k,v = x.split("=",1)
            d[k] = v
    return d

n_del = n_ins = n_ins_written = n_ins_skipped = 0
used = {}

with opener(vcf) as f, open(delbed,"w") as d, open(insfa,"w") as fa, open(insinfo,"w") as it:
    it.write("ID\tCHROM\tPOS\tSVLEN\n")

    for recno,line in enumerate(f,1):
        if line.startswith("#"): continue
        x = line.rstrip().split("\t")
        chrom, pos, vid, ref, alt = x[0], int(x[1]), x[2], x[3], x[4]
        inf = info_dict(x[7])
        typ = inf.get("SVTYPE","")
        svlen = abs(int(inf.get("SVLEN","0").split(",")[0])) if "SVLEN" in inf else 0

        key = vid if vid != "." else f"{typ}_{chrom}_{pos}"
        used[key] = used.get(key,0) + 1
        if used[key] > 1: key = f"{key}_{used[key]}"

        if typ == "DEL":
            n_del += 1
            end = int(inf["END"]) if "END" in inf else pos + svlen
            start = pos                    # excludes retained VCF anchor base
            length = end - start
            if length > 0:
                d.write(f"{chrom}\t{start}\t{end}\t{key}\t{length}\tDEL\n")

        elif typ == "INS":
            n_ins += 1

            # SVIM-asm normally stores explicit insertion sequence in ALT.
            if alt.startswith("<") or "[" in alt or "]" in alt or alt == ".":
                n_ins_skipped += 1
                continue

            if alt.startswith(ref) and len(alt) > len(ref):
                seq = alt[len(ref):]
            elif len(ref) == 1 and len(alt) > 1:
                seq = alt[1:]
            else:
                n_ins_skipped += 1
                continue

            if not seq:
                n_ins_skipped += 1
                continue

            n_ins_written += 1
            fa.write(f">{key}\n{seq}\n")
            it.write(f"{key}\t{chrom}\t{pos}\t{len(seq)}\n")

print(f"DEL records: {n_del}")
print(f"INS records: {n_ins}")
print(f"INS sequences extracted: {n_ins_written}")
print(f"INS sequences skipped: {n_ins_skipped}")
PY

# ============================================================
# 2. Convert Acra3RX RepeatMasker .out to BED
# ============================================================

python3 - "$RM_REF" "$REF_RM_BED" <<'PY'
import sys

inp, out = sys.argv[1:]

def broad(c):
    c = c.strip()
    u = c.upper()

    if u.startswith("LTR/COPIA"): return "LTR_Copia"
    if u.startswith("LTR/GYPSY") or u.startswith("LTR/DIRS"): return "LTR_Gypsy"
    if u.startswith("LTR/"): return "LTR_other"
    if u.startswith("LINE/") or u == "LINE": return "LINE"
    if u.startswith("SINE/") or u == "SINE": return "SINE"
    if u.startswith("DNA/") or u == "DNA": return "DNA_transposon"
    if u.startswith("RC/") or "HELITRON" in u: return "Rolling_circle"
    if "SIMPLE_REPEAT" in u: return "Simple_repeat"
    if "LOW_COMPLEXITY" in u: return "Low_complexity"
    if "SATELLITE" in u: return "Satellite"
    if u in ("UNKNOWN","UNCLASSIFIED") or u.startswith("UNKNOWN/"): return "Unclassified"
    return "Other"

n = 0
with open(inp) as f, open(out,"w") as o:
    for line in f:
        z = line.split()
        if len(z) < 11 or not z[0].isdigit(): continue

        chrom = z[4]
        start = int(z[5]) - 1
        end = int(z[6])
        repeat = z[9]
        rawclass = z[10]
        o.write(f"{chrom}\t{start}\t{end}\t{broad(rawclass)}\t{repeat}\t{rawclass}\n")
        n += 1

print("RepeatMasker annotations converted:", n)
PY

# ============================================================
# 3. Intersect supported DELs with Acra3RX repeats
# ============================================================

bedtools intersect -a "$DEL_BED" -b "$REF_RM_BED" -wo > "$DEL_INT"

# ============================================================
# 4. Run RepeatMasker on actual inserted sequences
# ============================================================

INS_RM_DIR="$OUTDIR/insertion_repeatmasker"
INS_RM_OUT="$INS_RM_DIR/$(basename "$INS_FA").out"

if [[ ! -s "$INS_RM_OUT" ]]; then
    echo "Running RepeatMasker on inserted sequences..."

    RepeatMasker \
        -pa "$THREADS" \
        -lib "$LIB" \
        -dir "$INS_RM_DIR" \
        "$INS_FA"
else
    echo "Existing insertion RepeatMasker output found; reusing it."
fi

# ============================================================
# 5. Summarize DEL and INS repeat associations
# ============================================================

python3 - "$DEL_BED" "$DEL_INT" "$INS_INFO" "$INS_RM_OUT" "$OUTDIR" <<'PY'
import sys
from collections import defaultdict, Counter

delbed, delint, insinfo, insrm, outdir = sys.argv[1:]

def merge_len(intervals):
    if not intervals: return 0
    intervals = sorted(intervals)
    total = 0
    s,e = intervals[0]
    for ns,ne in intervals[1:]:
        if ns <= e:
            e = max(e,ne)
        else:
            total += e-s
            s,e = ns,ne
    return total + e-s

def broad(c):
    u = c.strip().upper()
    if u.startswith("LTR/COPIA"): return "LTR_Copia"
    if u.startswith("LTR/GYPSY") or u.startswith("LTR/DIRS"): return "LTR_Gypsy"
    if u.startswith("LTR/"): return "LTR_other"
    if u.startswith("LINE/") or u == "LINE": return "LINE"
    if u.startswith("SINE/") or u == "SINE": return "SINE"
    if u.startswith("DNA/") or u == "DNA": return "DNA_transposon"
    if u.startswith("RC/") or "HELITRON" in u: return "Rolling_circle"
    if "SIMPLE_REPEAT" in u: return "Simple_repeat"
    if "LOW_COMPLEXITY" in u: return "Low_complexity"
    if "SATELLITE" in u: return "Satellite"
    if u in ("UNKNOWN","UNCLASSIFIED") or u.startswith("UNKNOWN/"): return "Unclassified"
    return "Other"

# ------------------------------------------------------------
# DEL
# ------------------------------------------------------------

dellen = {}
with open(delbed) as f:
    for line in f:
        z = line.rstrip().split("\t")
        dellen[z[3]] = int(z[4])

del_seg = defaultdict(lambda: defaultdict(list))
del_any = defaultdict(list)

with open(delint) as f:
    for line in f:
        z = line.rstrip().split("\t")
        aid = z[3]
        astart, aend = int(z[1]), int(z[2])
        bstart, bend = int(z[7]), int(z[8])
        cls = z[9]

        s = max(astart,bstart)
        e = min(aend,bend)
        if e > s:
            del_seg[aid][cls].append((s,e))
            del_any[aid].append((s,e))

del_event = []
del_class_events = defaultdict(int)
del_class_bp = defaultdict(int)
del_dominant = Counter()

for aid,length in dellen.items():
    perclass = {c:merge_len(v) for c,v in del_seg[aid].items()}
    anybp = merge_len(del_any[aid])

    for c,bp in perclass.items():
        if bp > 0:
            del_class_events[c] += 1
            del_class_bp[c] += bp

    dominant = max(perclass,key=perclass.get) if perclass else "Non_repeat"
    dominant_bp = perclass.get(dominant,0)
    del_dominant[dominant] += 1
    del_event.append((aid,length,anybp,100*anybp/length if length else 0,dominant,dominant_bp))

with open(f"{outdir}/deletion_repeat_by_event.tsv","w") as o:
    o.write("ID\tSVLEN\tREPEAT_BP\tREPEAT_PERCENT\tDOMINANT_CLASS\tDOMINANT_BP\n")
    for r in del_event: o.write("\t".join(map(str,r))+"\n")

with open(f"{outdir}/deletion_repeat_summary.tsv","w") as o:
    o.write("CLASS\tEVENTS_OVERLAPPING\tDOMINANT_EVENTS\tOVERLAP_BP\n")
    classes = sorted(set(del_class_events) | set(del_dominant))
    for c in classes:
        o.write(f"{c}\t{del_class_events.get(c,0)}\t{del_dominant.get(c,0)}\t{del_class_bp.get(c,0)}\n")

# ------------------------------------------------------------
# INS
# ------------------------------------------------------------

inslen = {}
with open(insinfo) as f:
    next(f)
    for line in f:
        z = line.rstrip().split("\t")
        inslen[z[0]] = int(z[3])

ins_seg = defaultdict(lambda: defaultdict(list))
ins_any = defaultdict(list)

with open(insrm) as f:
    for line in f:
        z = line.split()
        if len(z) < 11 or not z[0].isdigit(): continue

        seqid = z[4]
        start = int(z[5])-1
        end = int(z[6])
        cls = broad(z[10])

        ins_seg[seqid][cls].append((start,end))
        ins_any[seqid].append((start,end))

ins_event = []
ins_class_events = defaultdict(int)
ins_class_bp = defaultdict(int)
ins_dominant = Counter()

for aid,length in inslen.items():
    perclass = {c:merge_len(v) for c,v in ins_seg[aid].items()}
    anybp = merge_len(ins_any[aid])

    for c,bp in perclass.items():
        if bp > 0:
            ins_class_events[c] += 1
            ins_class_bp[c] += bp

    dominant = max(perclass,key=perclass.get) if perclass else "Non_repeat"
    dominant_bp = perclass.get(dominant,0)
    ins_dominant[dominant] += 1
    ins_event.append((aid,length,anybp,100*anybp/length if length else 0,dominant,dominant_bp))

with open(f"{outdir}/insertion_repeat_by_event.tsv","w") as o:
    o.write("ID\tSVLEN\tREPEAT_BP\tREPEAT_PERCENT\tDOMINANT_CLASS\tDOMINANT_BP\n")
    for r in ins_event: o.write("\t".join(map(str,r))+"\n")

with open(f"{outdir}/insertion_repeat_summary.tsv","w") as o:
    o.write("CLASS\tEVENTS_OVERLAPPING\tDOMINANT_EVENTS\tOVERLAP_BP\n")
    classes = sorted(set(ins_class_events) | set(ins_dominant))
    for c in classes:
        o.write(f"{c}\t{ins_class_events.get(c,0)}\t{ins_dominant.get(c,0)}\t{ins_class_bp.get(c,0)}\n")

# ------------------------------------------------------------
# Overall overview
# ------------------------------------------------------------

del_total = len(dellen)
del_rep = sum(1 for x in del_event if x[2] > 0)
del_bp = sum(dellen.values())
del_rep_bp = sum(merge_len(del_any[x]) for x in dellen)

ins_total = len(inslen)
ins_rep = sum(1 for x in ins_event if x[2] > 0)
ins_bp = sum(inslen.values())
ins_rep_bp = sum(merge_len(ins_any[x]) for x in inslen)

with open(f"{outdir}/sv_repeat_overview.tsv","w") as o:
    o.write("SVTYPE\tTOTAL_EVENTS\tREPEAT_ASSOCIATED\tEVENT_PERCENT\tTOTAL_BP\tREPEAT_BP\tBP_PERCENT\n")
    o.write(f"DEL\t{del_total}\t{del_rep}\t{100*del_rep/del_total:.2f}\t{del_bp}\t{del_rep_bp}\t{100*del_rep_bp/del_bp:.2f}\n")
    o.write(f"INS\t{ins_total}\t{ins_rep}\t{100*ins_rep/ins_total:.2f}\t{ins_bp}\t{ins_rep_bp}\t{100*ins_rep_bp/ins_bp:.2f}\n")

print("\n============================")
print("SV-repeat overview")
print("============================")
print(f"DEL: {del_rep}/{del_total} ({100*del_rep/del_total:.2f}%) overlap repeats")
print(f"     {100*del_rep_bp/del_bp:.2f}% of deleted sequence is repetitive")
print(f"INS: {ins_rep}/{ins_total} ({100*ins_rep/ins_total:.2f}%) contain repeats")
print(f"     {100*ins_rep_bp/ins_bp:.2f}% of inserted sequence is repetitive")
PY

echo
echo "============================================================"
echo "SV-repeat analysis complete"
echo "============================================================"
echo
cat "$OUTDIR/sv_repeat_overview.tsv"

echo
echo "--- Deletions ---"
cat "$OUTDIR/deletion_repeat_summary.tsv"

echo
echo "--- Insertions ---"
cat "$OUTDIR/insertion_repeat_summary.tsv"