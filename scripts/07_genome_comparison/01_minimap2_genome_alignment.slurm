#!/bin/bash
#SBATCH --job-name=acra_wga
#SBATCH --output=acra_wga_%j.out
#SBATCH --error=acra_wga_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

# ============================================================
# Whole-genome alignment + coverage + identity + dotplot
#
# Usage:
# sbatch acra_dotplot_align.slurm REF.fa QUERY.fa OUTDIR [MINLEN]
#
# REF   = Acra3RX reference genome
# QUERY = AcraUNRI genome
# MINLEN  = minimum PAF alignment block
#           default = 50000 bp
# ============================================================

REF="$1"
QUERY="$2"
OUTDIR="$3"
MINLEN="${4:-50000}"

THREADS="${SLURM_CPUS_PER_TASK:-32}"

MINIMAP2="${MINIMAP2:-minimap2}"
SAMTOOLS="${SAMTOOLS:-samtools}"

command -v "$MINIMAP2" >/dev/null 2>&1 || {
    echo "ERROR: minimap2 not found." >&2
    exit 1
}

command -v "$SAMTOOLS" >/dev/null 2>&1 || {
    echo "ERROR: samtools not found." >&2
    exit 1
}

mkdir -p "$OUTDIR"

RAW="$OUTDIR/acra_v4_vs_acra3rx.raw.paf"
PAF="$OUTDIR/acra_v4_vs_acra3rx.min${MINLEN}.paf"

SUMMARY="$OUTDIR/alignment_summary.min${MINLEN}.txt"
IDENTITY="$OUTDIR/alignment_identity.min${MINLEN}.tsv"
REFCOV="$OUTDIR/reference_scaffold_coverage.min${MINLEN}.tsv"
QUERYCOV="$OUTDIR/query_contig_coverage.min${MINLEN}.tsv"

REFBED="$OUTDIR/reference_exact_aligned.min${MINLEN}.bed"
QUERYBED="$OUTDIR/query_exact_aligned.min${MINLEN}.bed"

ORDER="$OUTDIR/query_dotplot_order.min${MINLEN}.tsv"
DOTPLOT="$OUTDIR/acra_v4_vs_acra3rx.min${MINLEN}.dotplot.svg"

echo "========================================================"
echo "A. crassicarpa whole-genome comparison"
echo "========================================================"
echo "Reference : $REF"
echo "Query     : $QUERY"
echo "Min block : $MINLEN bp"
echo

[ -f "${REF}.fai" ] || "$SAMTOOLS" faidx "$REF"
[ -f "${QUERY}.fai" ] || "$SAMTOOLS" faidx "$QUERY"

# ============================================================
# 1. Whole-genome alignment
# ============================================================

echo "Running minimap2..."

"$MINIMAP2" \
    -x asm5 \
    -c \
    --eqx \
    --secondary=no \
    -t "${THREADS}" \
    "$REF" "$QUERY" \
    > "$RAW"

echo "Raw alignments: $(wc -l < "$RAW")"

# ============================================================
# 2. Filter
# ============================================================

awk -v m="$MINLEN" '$11 >= m' "$RAW" > "$PAF"

echo "Alignments >= ${MINLEN} bp: $(wc -l < "$PAF")"

# ============================================================
# 3. Statistics + exact coverage + dotplot
# ============================================================

python3 - \
"$PAF" \
"${REF}.fai" \
"${QUERY}.fai" \
"$SUMMARY" \
"$IDENTITY" \
"$REFCOV" \
"$QUERYCOV" \
"$REFBED" \
"$QUERYBED" \
"$ORDER" \
"$DOTPLOT" \
"$MINLEN" <<'PY'

import sys, re, statistics
from collections import defaultdict

(
    paf_file,
    ref_fai,
    query_fai,
    summary_file,
    identity_file,
    refcov_file,
    querycov_file,
    refbed_file,
    querybed_file,
    order_file,
    svg_file,
    minlen
) = sys.argv[1:]

minlen = int(minlen)

# ============================================================
# FASTA indices
# ============================================================

def read_fai(fn):
    order = []
    lengths = {}

    with open(fn) as f:
        for line in f:
            x = line.split()
            order.append(x[0])
            lengths[x[0]] = int(x[1])

    return order, lengths

ref_order, ref_len = read_fai(ref_fai)
query_order_original, query_len = read_fai(query_fai)

ref_index = {x:i for i,x in enumerate(ref_order)}

# ============================================================
# Helpers
# ============================================================

cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')

def get_cigar(fields):

    for z in fields[12:]:
        if z.startswith("cg:Z:"):
            return z[5:]

    return None


def merge_intervals(intervals):

    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for s,e in intervals[1:]:

        if s <= merged[-1][1]:
            merged[-1][1] = max(
                merged[-1][1], e
            )

        else:
            merged.append([s,e])

    return [(s,e) for s,e in merged]

# ============================================================
# Parse alignments
# ============================================================

alignments = []

ref_intervals = defaultdict(list)
query_intervals = defaultdict(list)

anchor_bp = defaultdict(lambda: defaultdict(int))
anchor_pos_num = defaultdict(lambda: defaultdict(float))
anchor_pos_den = defaultdict(lambda: defaultdict(float))
strand_bp = defaultdict(
    lambda: defaultdict(
        lambda: defaultdict(int)
    )
)

with open(paf_file) as f:

    for line in f:

        x = line.rstrip().split("\t")

        qname = x[0]
        qlen = int(x[1])
        qstart = int(x[2])
        qend = int(x[3])

        strand = x[4]

        tname = x[5]
        tlen = int(x[6])
        tstart = int(x[7])
        tend = int(x[8])

        paf_matches = int(x[9])
        paf_alnlen = int(x[10])
        mapq = int(x[11])

        cigar = get_cigar(x)

        if cigar is None:
            raise RuntimeError(
                "No cg:Z tag found. "
                "minimap2 must be run with -c."
            )

        eq_bp = 0
        x_bp = 0
        ins_bp = 0
        del_bp = 0

        rpos = tstart

        if strand == "+":
            qpos = qstart
        else:
            qpos = qend

        # ----------------------------------------------------
        # CIGAR
        # ----------------------------------------------------

        for n,op in cigar_re.findall(cigar):

            n = int(n)

            if op == "=":

                eq_bp += n

                ref_intervals[tname].append(
                    (rpos, rpos+n)
                )

                if strand == "+":
                    query_intervals[qname].append(
                        (qpos, qpos+n)
                    )
                    qpos += n
                else:
                    query_intervals[qname].append(
                        (qpos-n, qpos)
                    )
                    qpos -= n

                rpos += n

            elif op == "X":

                x_bp += n

                ref_intervals[tname].append(
                    (rpos, rpos+n)
                )

                if strand == "+":
                    query_intervals[qname].append(
                        (qpos, qpos+n)
                    )
                    qpos += n
                else:
                    query_intervals[qname].append(
                        (qpos-n, qpos)
                    )
                    qpos -= n

                rpos += n

            elif op == "M":

                # fallback if M occurs
                ref_intervals[tname].append(
                    (rpos, rpos+n)
                )

                if strand == "+":
                    query_intervals[qname].append(
                        (qpos, qpos+n)
                    )
                    qpos += n
                else:
                    query_intervals[qname].append(
                        (qpos-n, qpos)
                    )
                    qpos -= n

                rpos += n

            elif op == "I":

                ins_bp += n

                if strand == "+":
                    qpos += n
                else:
                    qpos -= n

            elif op in ("D","N"):

                del_bp += n
                rpos += n

            elif op == "S":

                if strand == "+":
                    qpos += n
                else:
                    qpos -= n

            elif op in ("H","P"):
                pass

        aligned_pair_bp = eq_bp + x_bp

        nucleotide_identity = (
            100.0 * eq_bp / aligned_pair_bp
            if aligned_pair_bp > 0 else 0
        )

        gap_denom = (
            eq_bp +
            x_bp +
            ins_bp +
            del_bp
        )

        gap_identity = (
            100.0 * eq_bp / gap_denom
            if gap_denom > 0 else 0
        )

        a = {
            "qname":qname,
            "qlen":qlen,
            "qstart":qstart,
            "qend":qend,
            "strand":strand,

            "tname":tname,
            "tlen":tlen,
            "tstart":tstart,
            "tend":tend,

            "paf_alnlen":paf_alnlen,
            "mapq":mapq,

            "eq_bp":eq_bp,
            "x_bp":x_bp,
            "ins_bp":ins_bp,
            "del_bp":del_bp,

            "aligned_pair_bp":aligned_pair_bp,

            "nuc_identity":nucleotide_identity,
            "gap_identity":gap_identity
        }

        alignments.append(a)

        # ----------------------------------------------------
        # placement of query contigs
        # ----------------------------------------------------

        anchor_bp[qname][tname] += aligned_pair_bp

        midpoint = (tstart+tend)/2

        anchor_pos_num[qname][tname] += (
            midpoint * aligned_pair_bp
        )

        anchor_pos_den[qname][tname] += (
            aligned_pair_bp
        )

        strand_bp[qname][tname][strand] += (
            aligned_pair_bp
        )

# ============================================================
# Merge coverage
# ============================================================

ref_merged = {
    c: merge_intervals(v)
    for c,v in ref_intervals.items()
}

query_merged = {
    c: merge_intervals(v)
    for c,v in query_intervals.items()
}

# ============================================================
# BED
# ============================================================

with open(refbed_file,"w") as o:

    for chrom in ref_order:

        for s,e in ref_merged.get(chrom,[]):
            o.write(
                f"{chrom}\t{s}\t{e}\n"
            )


with open(querybed_file,"w") as o:

    for chrom in query_order_original:

        for s,e in query_merged.get(chrom,[]):
            o.write(
                f"{chrom}\t{s}\t{e}\n"
            )

# ============================================================
# Exact genome coverage
# ============================================================

ref_total = sum(ref_len.values())
query_total = sum(query_len.values())

ref_cov_bp = sum(
    e-s
    for c in ref_merged
    for s,e in ref_merged[c]
)

query_cov_bp = sum(
    e-s
    for c in query_merged
    for s,e in query_merged[c]
)

ref_cov_pct = (
    100*ref_cov_bp/ref_total
)

query_cov_pct = (
    100*query_cov_bp/query_total
)

# ============================================================
# Reference scaffold coverage
# ============================================================

with open(refcov_file,"w") as o:

    o.write(
        "SCAFFOLD\tLENGTH_BP\t"
        "ALIGNED_BP\tCOVERAGE_PCT\n"
    )

    for chrom in ref_order:

        cov = sum(
            e-s
            for s,e
            in ref_merged.get(chrom,[])
        )

        pct = (
            100*cov/ref_len[chrom]
        )

        o.write(
            f"{chrom}\t"
            f"{ref_len[chrom]}\t"
            f"{cov}\t"
            f"{pct:.3f}\n"
        )

# ============================================================
# Query contig coverage
# ============================================================

with open(querycov_file,"w") as o:

    o.write(
        "CONTIG\tLENGTH_BP\t"
        "ALIGNED_BP\tCOVERAGE_PCT\n"
    )

    for chrom in query_order_original:

        cov = sum(
            e-s
            for s,e
            in query_merged.get(chrom,[])
        )

        pct = (
            100*cov/query_len[chrom]
        )

        o.write(
            f"{chrom}\t"
            f"{query_len[chrom]}\t"
            f"{cov}\t"
            f"{pct:.3f}\n"
        )

# ============================================================
# Identity output
# ============================================================

with open(identity_file,"w") as o:

    o.write(
        "QUERY\tREFERENCE\t"
        "QUERY_START\tQUERY_END\t"
        "REF_START\tREF_END\t"
        "STRAND\t"
        "MATCH_BP\tMISMATCH_BP\t"
        "INS_BP\tDEL_BP\t"
        "NUCLEOTIDE_IDENTITY_PCT\t"
        "GAP_INCLUSIVE_IDENTITY_PCT\t"
        "MAPQ\n"
    )

    for a in alignments:

        o.write(
            f"{a['qname']}\t"
            f"{a['tname']}\t"
            f"{a['qstart']}\t"
            f"{a['qend']}\t"
            f"{a['tstart']}\t"
            f"{a['tend']}\t"
            f"{a['strand']}\t"
            f"{a['eq_bp']}\t"
            f"{a['x_bp']}\t"
            f"{a['ins_bp']}\t"
            f"{a['del_bp']}\t"
            f"{a['nuc_identity']:.4f}\t"
            f"{a['gap_identity']:.4f}\t"
            f"{a['mapq']}\n"
        )

# ============================================================
# Correct identity statistics
# ============================================================

TOTAL_EQ = sum(
    a["eq_bp"] for a in alignments
)

TOTAL_X = sum(
    a["x_bp"] for a in alignments
)

TOTAL_INS = sum(
    a["ins_bp"] for a in alignments
)

TOTAL_DEL = sum(
    a["del_bp"] for a in alignments
)

TOTAL_PAIR = TOTAL_EQ + TOTAL_X

weighted_nuc_identity = (
    100*TOTAL_EQ/TOTAL_PAIR
)

weighted_gap_identity = (
    100*TOTAL_EQ /
    (
        TOTAL_EQ +
        TOTAL_X +
        TOTAL_INS +
        TOTAL_DEL
    )
)

nuc_ids = [
    a["nuc_identity"]
    for a in alignments
    if a["aligned_pair_bp"] > 0
]

mean_nuc_identity = statistics.mean(
    nuc_ids
)

median_nuc_identity = statistics.median(
    nuc_ids
)

n95 = sum(
    x >= 95 for x in nuc_ids
)

n99 = sum(
    x >= 99 for x in nuc_ids
)

# ============================================================
# Query contig ordering
# ============================================================

query_anchor = {}

for qname in anchor_bp:

    best_ref = max(
        anchor_bp[qname],
        key=lambda r:
        anchor_bp[qname][r]
    )

    pos = (
        anchor_pos_num[qname][best_ref] /
        anchor_pos_den[qname][best_ref]
    )

    plus = (
        strand_bp[qname][best_ref]["+"]
    )

    minus = (
        strand_bp[qname][best_ref]["-"]
    )

    orientation = (
        "+" if plus >= minus else "-"
    )

    query_anchor[qname] = (
        best_ref,
        pos,
        orientation,
        anchor_bp[qname][best_ref]
    )

ordered_query = sorted(
    query_anchor,
    key=lambda q: (
        ref_index.get(
            query_anchor[q][0],
            999999
        ),
        query_anchor[q][1]
    )
)

unaligned = [
    q for q in query_order_original
    if q not in query_anchor
]

ordered_query.extend(unaligned)

with open(order_file,"w") as o:

    o.write(
        "CONTIG\tLENGTH_BP\t"
        "ANCHOR_REFERENCE\t"
        "ANCHOR_POSITION\t"
        "ORIENTATION\t"
        "ANCHOR_ALIGNED_BP\n"
    )

    for q in ordered_query:

        if q in query_anchor:

            r,p,ori,bp = (
                query_anchor[q]
            )

            o.write(
                f"{q}\t"
                f"{query_len[q]}\t"
                f"{r}\t"
                f"{p:.1f}\t"
                f"{ori}\t"
                f"{bp}\n"
            )

        else:

            o.write(
                f"{q}\t"
                f"{query_len[q]}\t"
                f"NA\tNA\tNA\t0\n"
            )

# ============================================================
# Cumulative coordinates
# ============================================================

ref_offset = {}
x = 0

for chrom in ref_order:

    ref_offset[chrom] = x
    x += ref_len[chrom]

ref_plot_total = x

query_offset = {}
y = 0

for q in ordered_query:

    query_offset[q] = y
    y += query_len[q]

query_plot_total = y

# ============================================================
# SVG dotplot
# ============================================================

W = 1600
H = 1250

LEFT = 125
RIGHT = 55
TOP = 75
BOTTOM = 120

PW = W-LEFT-RIGHT
PH = H-TOP-BOTTOM

def sx(bp):
    return (
        LEFT +
        PW*bp/ref_plot_total
    )

def sy(bp):
    return (
        TOP +
        PH -
        PH*bp/query_plot_total
    )

with open(svg_file,"w") as o:

    o.write(
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}">\n'
    )

    o.write(
        '<rect width="100%" height="100%" '
        'fill="white"/>\n'
    )

    # Title

    o.write(
        '<text x="800" y="35" '
        'text-anchor="middle" '
        'font-family="Arial,sans-serif" '
        'font-size="26" '
        'font-weight="bold">'
        'Whole-genome alignment of '
        'Acacia crassicarpa assemblies'
        '</text>\n'
    )

    o.write(
        f'<text x="800" y="60" '
        f'text-anchor="middle" '
        f'font-family="Arial,sans-serif" '
        f'font-size="14" '
        f'fill="#555555">'
        f'Alignment blocks ≥ '
        f'{minlen/1000:g} kb'
        f'</text>\n'
    )

    # Border

    o.write(
        f'<rect x="{LEFT}" y="{TOP}" '
        f'width="{PW}" height="{PH}" '
        f'fill="none" '
        f'stroke="black" '
        f'stroke-width="1.2"/>\n'
    )

    # --------------------------------------------------------
    # Reference boundaries
    #
    # Draw only scaffolds >=1 Mb
    # --------------------------------------------------------

    cumulative = 0

    for chrom in ref_order:

        L = ref_len[chrom]

        if (
            cumulative > 0
            and L >= 1_000_000
        ):

            xx = sx(cumulative)

            o.write(
                f'<line '
                f'x1="{xx:.2f}" '
                f'y1="{TOP}" '
                f'x2="{xx:.2f}" '
                f'y2="{TOP+PH}" '
                f'stroke="#dddddd" '
                f'stroke-width="0.7"/>\n'
            )

        cumulative += L

    # --------------------------------------------------------
    # Alignment segments
    # --------------------------------------------------------

    for a in alignments:

        q = a["qname"]

        if q not in query_offset:
            continue

        qlen = query_len[q]

        orientation = (
            query_anchor[q][2]
            if q in query_anchor
            else "+"
        )

        if a["strand"] == "+":
            qa = a["qstart"]
            qb = a["qend"]

        else:
            qa = a["qend"]
            qb = a["qstart"]

        if orientation == "-":
            qa = qlen-qa
            qb = qlen-qb

        same_orientation = (
            a["strand"] == orientation
        )

        colour = (
            "#234E70"
            if same_orientation
            else "#B23A48"
        )

        x1 = sx(
            ref_offset[a["tname"]]
            + a["tstart"]
        )

        x2 = sx(
            ref_offset[a["tname"]]
            + a["tend"]
        )

        y1 = sy(
            query_offset[q] + qa
        )

        y2 = sy(
            query_offset[q] + qb
        )

        o.write(
            f'<line '
            f'x1="{x1:.2f}" '
            f'y1="{y1:.2f}" '
            f'x2="{x2:.2f}" '
            f'y2="{y2:.2f}" '
            f'stroke="{colour}" '
            f'stroke-width="0.9" '
            f'stroke-opacity="0.75"/>\n'
        )

    # --------------------------------------------------------
    # Axis labels
    # --------------------------------------------------------

    o.write(
        f'<text '
        f'x="{LEFT+PW/2}" '
        f'y="{H-28}" '
        f'text-anchor="middle" '
        f'font-family="Arial,sans-serif" '
        f'font-size="20">'
        f'Acra3RX reference '
        f'(cumulative Mb)'
        f'</text>\n'
    )

    o.write(
        f'<text '
        f'x="27" '
        f'y="{TOP+PH/2}" '
        f'text-anchor="middle" '
        f'font-family="Arial,sans-serif" '
        f'font-size="20" '
        f'transform="rotate(-90 27 '
        f'{TOP+PH/2})">'
        f'v4 assembly '
        f'(ordered contigs, cumulative Mb)'
        f'</text>\n'
    )

    # X ticks every 100 Mb

    tick = 0

    while tick <= ref_plot_total:

        xx = sx(tick)

        o.write(
            f'<line '
            f'x1="{xx:.2f}" '
            f'y1="{TOP+PH}" '
            f'x2="{xx:.2f}" '
            f'y2="{TOP+PH+7}" '
            f'stroke="black"/>\n'
        )

        o.write(
            f'<text '
            f'x="{xx:.2f}" '
            f'y="{TOP+PH+27}" '
            f'text-anchor="middle" '
            f'font-family="Arial,sans-serif" '
            f'font-size="13">'
            f'{tick/1e6:.0f}'
            f'</text>\n'
        )

        tick += 100_000_000

    # Y ticks every 100 Mb

    tick = 0

    while tick <= query_plot_total:

        yy = sy(tick)

        o.write(
            f'<line '
            f'x1="{LEFT-7}" '
            f'y1="{yy:.2f}" '
            f'x2="{LEFT}" '
            f'y2="{yy:.2f}" '
            f'stroke="black"/>\n'
        )

        o.write(
            f'<text '
            f'x="{LEFT-12}" '
            f'y="{yy+4:.2f}" '
            f'text-anchor="end" '
            f'font-family="Arial,sans-serif" '
            f'font-size="13">'
            f'{tick/1e6:.0f}'
            f'</text>\n'
        )

        tick += 100_000_000

    # Legend

    lx = W-320
    ly = 100

    o.write(
        f'<line '
        f'x1="{lx}" y1="{ly}" '
        f'x2="{lx+45}" y2="{ly}" '
        f'stroke="#234E70" '
        f'stroke-width="4"/>\n'
    )

    o.write(
        f'<text '
        f'x="{lx+58}" '
        f'y="{ly+5}" '
        f'font-family="Arial,sans-serif" '
        f'font-size="14">'
        f'Same orientation'
        f'</text>\n'
    )

    o.write(
        f'<line '
        f'x1="{lx}" y1="{ly+25}" '
        f'x2="{lx+45}" y2="{ly+25}" '
        f'stroke="#B23A48" '
        f'stroke-width="4"/>\n'
    )

    o.write(
        f'<text '
        f'x="{lx+58}" '
        f'y="{ly+30}" '
        f'font-family="Arial,sans-serif" '
        f'font-size="14">'
        f'Opposite orientation'
        f'</text>\n'
    )

    o.write('</svg>\n')

# ============================================================
# Summary
# ============================================================

N = len(alignments)

with open(summary_file,"w") as o:

    o.write(
        "==================================================\n"
    )

    o.write(
        "Acacia crassicarpa whole-genome alignment summary\n"
    )

    o.write(
        "==================================================\n\n"
    )

    o.write(
        f"Minimum alignment block: "
        f"{minlen} bp\n\n"
    )

    o.write("ASSEMBLY SIZE\n")

    o.write(
        f"Reference total bp: "
        f"{ref_total}\n"
    )

    o.write(
        f"Query total bp:     "
        f"{query_total}\n\n"
    )

    o.write("FILTERED ALIGNMENTS\n")

    o.write(
        f"Number of alignment blocks: "
        f"{N}\n\n"
    )

    o.write(
        "EXACT ALIGNED-BASE COVERAGE\n"
    )

    o.write(
        f"Reference exact aligned bp: "
        f"{ref_cov_bp}\n"
    )

    o.write(
        f"Reference exact coverage:   "
        f"{ref_cov_pct:.3f}%\n"
    )

    o.write(
        f"Query exact aligned bp:     "
        f"{query_cov_bp}\n"
    )

    o.write(
        f"Query exact coverage:       "
        f"{query_cov_pct:.3f}%\n\n"
    )

    o.write("ALIGNMENT IDENTITY\n")

    o.write(
        f"Aligned-base nucleotide identity: "
        f"{weighted_nuc_identity:.4f}%\n"
    )

    o.write(
        f"Gap-inclusive identity:            "
        f"{weighted_gap_identity:.4f}%\n"
    )

    o.write(
        f"Mean block nucleotide identity:    "
        f"{mean_nuc_identity:.4f}%\n"
    )

    o.write(
        f"Median block nucleotide identity:  "
        f"{median_nuc_identity:.4f}%\n"
    )

    o.write(
        f"Blocks >=95% nucleotide identity: "
        f"{n95} "
        f"({100*n95/N:.2f}%)\n"
    )

    o.write(
        f"Blocks >=99% nucleotide identity: "
        f"{n99} "
        f"({100*n99/N:.2f}%)\n\n"
    )

    o.write(
        "ALIGNMENT LENGTH DISTRIBUTION\n"
    )

    for threshold in (
        10000,
        50000,
        100000,
        500000,
        1000000
    ):

        n = sum(
            a["paf_alnlen"] >= threshold
            for a in alignments
        )

        o.write(
            f">={threshold/1000:g} kb: "
            f"{n}\n"
        )

    o.write("\nDOTPLOT\n")

    o.write(
        f"Aligned query contigs: "
        f"{len(query_anchor)}\n"
    )

    o.write(
        f"Unaligned query contigs: "
        f"{len(unaligned)}\n"
    )

    o.write(
        "Query contigs were ordered and oriented "
        "according to their dominant Acra3RX alignment.\n"
    )

print()
print(open(summary_file).read())

print("Dotplot:", svg_file)

PY

echo
echo "========================================================"
echo "Completed"
echo "========================================================"
echo "Summary:   $SUMMARY"
echo "Dotplot:   $DOTPLOT"
echo "Identity:  $IDENTITY"
echo "Ref cov:   $REFCOV"
echo "Query cov: $QUERYCOV"
echo