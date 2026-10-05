#!/usr/bin/env python3

import argparse, gzip, os, random, statistics
from collections import defaultdict
from bisect import bisect_left

p = argparse.ArgumentParser()
p.add_argument("--vcf", required=True)
p.add_argument("--gtf", required=True)
p.add_argument("--fai", required=True)
p.add_argument("--outdir", required=True)
p.add_argument("--nperm", type=int, default=1000)
p.add_argument("--seed", type=int, default=42)
a = p.parse_args()

os.makedirs(a.outdir, exist_ok=True)
random.seed(a.seed)

def opener(x):
    return gzip.open(x, "rt") if x.endswith(".gz") else open(x)

# ============================================================
# Reference scaffold lengths
# ============================================================

scaflen = {}
with open(a.fai) as f:
    for line in f:
        x = line.split()
        scaflen[x[0]] = int(x[1])

# ============================================================
# Read and merge GTF gene/CDS intervals
# GTF = 1-based inclusive -> BED-like 0-based half-open
# ============================================================

raw_gene = defaultdict(list)
raw_cds = defaultdict(list)

with open(a.gtf) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        x = line.rstrip().split(None, 8)
        if len(x) < 9:
            continue

        chrom, feature = x[0], x[2]
        start, end = int(x[3]) - 1, int(x[4])

        if feature == "gene":
            raw_gene[chrom].append((start, end))
        elif feature == "CDS":
            raw_cds[chrom].append((start, end))

def merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged = [list(intervals[0])]

    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])

    return [(s, e) for s, e in merged]

def make_index(raw):
    idx = {}

    for chrom, iv in raw.items():
        m = merge_intervals(iv)
        idx[chrom] = (
            [x[0] for x in m],
            [x[1] for x in m]
        )

    return idx

gene_idx = make_index(raw_gene)
cds_idx = make_index(raw_cds)

def overlaps(index, chrom, start, end):
    if chrom not in index:
        return False

    starts, ends = index[chrom]

    # rightmost interval starting before event end
    i = bisect_left(starts, end) - 1

    return i >= 0 and ends[i] > start

# ============================================================
# Read supported SVs
# DEL interval = sequence deleted after VCF anchor
# INS = breakpoint represented by the base immediately right
#       of the VCF anchor, consistent with previous analysis
# ============================================================

DEL = defaultdict(list)
INS = defaultdict(list)

observed = {
    ("DEL", "gene"): 0,
    ("DEL", "CDS"): 0,
    ("INS", "gene"): 0,
    ("INS", "CDS"): 0
}

def parse_info(s):
    d = {}
    for z in s.split(";"):
        if "=" in z:
            k, v = z.split("=", 1)
            d[k] = v
    return d

with opener(a.vcf) as f:
    for line in f:
        if line.startswith("#"):
            continue

        x = line.rstrip().split("\t")
        chrom = x[0]
        pos = int(x[1])
        info = parse_info(x[7])
        typ = info.get("SVTYPE", "")

        if typ == "DEL":
            if "END" in info:
                end = int(info["END"])
            else:
                svlen = abs(int(info["SVLEN"].split(",")[0]))
                end = pos + svlen

            start = pos
            length = end - start

            if length <= 0:
                continue

            DEL[chrom].append((start, length))

            if overlaps(gene_idx, chrom, start, end):
                observed[("DEL", "gene")] += 1

            if overlaps(cds_idx, chrom, start, end):
                observed[("DEL", "CDS")] += 1

        elif typ == "INS":
            start = pos
            end = pos + 1

            INS[chrom].append(start)

            if overlaps(gene_idx, chrom, start, end):
                observed[("INS", "gene")] += 1

            if overlaps(cds_idx, chrom, start, end):
                observed[("INS", "CDS")] += 1

n_del = sum(len(v) for v in DEL.values())
n_ins = sum(len(v) for v in INS.values())
n_all = n_del + n_ins

print("DEL:", n_del)
print("INS:", n_ins)
print("Total:", n_all)
print()
print("Observed gene-body overlap:")
print("DEL:", observed[("DEL", "gene")])
print("INS:", observed[("INS", "gene")])
print("ALL:", observed[("DEL", "gene")] + observed[("INS", "gene")])

# ============================================================
# Permutations
# ============================================================

null = {
    ("DEL", "gene"): [],
    ("DEL", "CDS"): [],
    ("INS", "gene"): [],
    ("INS", "CDS"): [],
    ("ALL", "gene"): [],
    ("ALL", "CDS"): []
}

null_file = os.path.join(a.outdir, "sv_gene_depletion_null.tsv")

with open(null_file, "w") as out:
    out.write("PERM\tDEL_GENE\tDEL_CDS\tINS_GENE\tINS_CDS\tALL_GENE\tALL_CDS\n")

    for perm in range(1, a.nperm + 1):

        dg = dc = ig = ic = 0

        # -------------------------
        # DEL: same scaffold + size
        # -------------------------
        for chrom, events in DEL.items():
            L = scaflen[chrom]

            for _, length in events:
                max_start = L - length

                if max_start < 0:
                    continue

                start = random.randint(0, max_start)
                end = start + length

                g = overlaps(gene_idx, chrom, start, end)

                if g:
                    dg += 1

                    # CDS is a subset of gene space, so only test if gene hit
                    if overlaps(cds_idx, chrom, start, end):
                        dc += 1

        # -------------------------
        # INS: same scaffold
        # -------------------------
        for chrom, events in INS.items():
            L = scaflen[chrom]

            for _ in events:
                # insertion breakpoint after a reference base
                start = random.randint(1, L)
                end = start + 1

                g = overlaps(gene_idx, chrom, start, end)

                if g:
                    ig += 1

                    if overlaps(cds_idx, chrom, start, end):
                        ic += 1

        ag = dg + ig
        ac = dc + ic

        null[("DEL", "gene")].append(dg)
        null[("DEL", "CDS")].append(dc)
        null[("INS", "gene")].append(ig)
        null[("INS", "CDS")].append(ic)
        null[("ALL", "gene")].append(ag)
        null[("ALL", "CDS")].append(ac)

        out.write(f"{perm}\t{dg}\t{dc}\t{ig}\t{ic}\t{ag}\t{ac}\n")

        if perm % 100 == 0:
            print(f"Completed {perm}/{a.nperm} permutations")

# ============================================================
# Statistics
# ============================================================

observed[("ALL", "gene")] = observed[("DEL", "gene")] + observed[("INS", "gene")]
observed[("ALL", "CDS")] = observed[("DEL", "CDS")] + observed[("INS", "CDS")]

totals = {
    "DEL": n_del,
    "INS": n_ins,
    "ALL": n_all
}

def percentile(x, q):
    x = sorted(x)

    if not x:
        return float("nan")

    k = (len(x)-1) * q
    lo = int(k)
    hi = min(lo+1, len(x)-1)
    frac = k-lo

    return x[lo]*(1-frac) + x[hi]*frac

rows = []

for typ in ("DEL", "INS", "ALL"):
    for feature in ("gene", "CDS"):

        obs = observed[(typ, feature)]
        vals = null[(typ, feature)]
        total = totals[typ]

        mean_null = statistics.mean(vals)
        sd_null = statistics.pstdev(vals)

        obs_pct = 100 * obs / total
        null_pct = 100 * mean_null / total

        low = 100 * percentile(vals, 0.025) / total
        high = 100 * percentile(vals, 0.975) / total

        fold = obs_pct / null_pct if null_pct else float("nan")

        # one-sided empirical depletion P
        pdep = (1 + sum(v <= obs for v in vals)) / (a.nperm + 1)

        rows.append({
            "type": typ,
            "feature": feature,
            "total": total,
            "obs": obs,
            "obs_pct": obs_pct,
            "null_mean": mean_null,
            "null_pct": null_pct,
            "ci_low": low,
            "ci_high": high,
            "fold": fold,
            "p": pdep
        })

# Benjamini-Hochberg FDR
m = len(rows)
order = sorted(range(m), key=lambda i: rows[i]["p"])

qvals = [0] * m
running = 1.0

for rank in range(m, 0, -1):
    i = order[rank-1]
    q = rows[i]["p"] * m / rank
    running = min(running, q)
    qvals[i] = min(running, 1.0)

summary = os.path.join(a.outdir, "sv_gene_depletion_summary.tsv")

with open(summary, "w") as out:
    out.write(
        "SVTYPE\tFEATURE\tTOTAL_EVENTS\tOBS_EVENTS\tOBS_PCT\t"
        "NULL_MEAN_EVENTS\tNULL_MEAN_PCT\tNULL_CI95_LOW\tNULL_CI95_HIGH\t"
        "OBS_NULL_FOLD\tP_DEPLETION\tFDR_DEPLETION\n"
    )

    for i, r in enumerate(rows):
        out.write(
            f"{r['type']}\t{r['feature']}\t{r['total']}\t{r['obs']}\t"
            f"{r['obs_pct']:.3f}\t{r['null_mean']:.3f}\t{r['null_pct']:.3f}\t"
            f"{r['ci_low']:.3f}\t{r['ci_high']:.3f}\t"
            f"{r['fold']:.3f}\t{r['p']:.6g}\t{qvals[i]:.6g}\n"
        )

print()
print("==========================================")
print("RESULTS")
print("==========================================")

with open(summary) as f:
    print(f.read())