#!/usr/bin/env python3

import os, argparse
import numpy as np
from collections import defaultdict

p = argparse.ArgumentParser()
p.add_argument("--deletions", required=True, help="supported_deletions.bed")
p.add_argument("--repeats", required=True, help="acra3rx_repeats.bed")
p.add_argument("--fai", required=True, help="Acra3RX FASTA .fai")
p.add_argument("--outdir", default="output")
p.add_argument("--nperm", type=int, default=1000)
p.add_argument("--seed", type=int, default=42)
a = p.parse_args()

os.makedirs(a.outdir, exist_ok=True)
rng = np.random.default_rng(a.seed)

# ---------- reference lengths ----------
genome = {}
with open(a.fai) as f:
    for line in f:
        x = line.split()
        genome[x[0]] = int(x[1])

# ---------- deletion intervals ----------
dels = defaultdict(list)
with open(a.deletions) as f:
    for line in f:
        x = line.rstrip().split("\t")
        chrom, start, end = x[0], int(x[1]), int(x[2])
        if chrom in genome and end > start:
            dels[chrom].append((start, end))

dels = {c: np.array(v, dtype=np.int64) for c, v in dels.items()}
total_events = sum(len(v) for v in dels.values())
total_bp = sum(np.sum(v[:,1] - v[:,0]) for v in dels.values())

# ---------- repeat intervals by class ----------
raw = defaultdict(lambda: defaultdict(list))
allrep = defaultdict(list)

with open(a.repeats) as f:
    for line in f:
        x = line.rstrip().split("\t")
        chrom, start, end, cls = x[0], int(x[1]), int(x[2]), x[3]
        if chrom in genome and end > start:
            raw[cls][chrom].append((start, end))
            allrep[chrom].append((start, end))

def merge_intervals(v):
    if not v: return np.empty((0,2), dtype=np.int64)
    v = sorted(v)
    out = [list(v[0])]
    for s,e in v[1:]:
        if s <= out[-1][1]: out[-1][1] = max(out[-1][1], e)
        else: out.append([s,e])
    return np.array(out, dtype=np.int64)

raw["Any_repeat"] = allrep
classes = ["Any_repeat"] + sorted(c for c in raw if c != "Any_repeat")

# Precompute merged repeat intervals and cumulative coverage.
models = {}
for cls in classes:
    models[cls] = {}
    for chrom in genome:
        iv = merge_intervals(raw[cls].get(chrom, []))
        if len(iv):
            starts, ends = iv[:,0], iv[:,1]
            prefix = np.concatenate(([0], np.cumsum(ends-starts)))
            models[cls][chrom] = (starts, ends, prefix)

def covered_to(x, model):
    starts, ends, prefix = model
    idx = np.searchsorted(starts, x, side="left") - 1
    out = np.zeros(len(x), dtype=np.int64)
    good = idx >= 0
    j = idx[good]
    out[good] = prefix[j] + np.minimum(np.maximum(x[good]-starts[j], 0), ends[j]-starts[j])
    return out

def overlap_bp(starts, ends, model):
    if model is None: return np.zeros(len(starts), dtype=np.int64)
    return covered_to(ends, model) - covered_to(starts, model)

# ---------- observed ----------
obs_bp = np.zeros(len(classes), dtype=np.int64)
obs_events = np.zeros(len(classes), dtype=np.int64)

for chrom, q in dels.items():
    for i, cls in enumerate(classes):
        ov = overlap_bp(q[:,0], q[:,1], models[cls].get(chrom))
        obs_bp[i] += ov.sum()
        obs_events[i] += np.sum(ov > 0)

print(f"Deletion events: {total_events}")
print(f"Deletion bp:     {total_bp}")
print(f"Permutations:    {a.nperm}")

# ---------- permutations ----------
null_bp = np.zeros((a.nperm, len(classes)), dtype=np.int64)
null_events = np.zeros((a.nperm, len(classes)), dtype=np.int64)

for perm in range(a.nperm):
    for chrom, q in dels.items():
        lengths = q[:,1] - q[:,0]
        maxstart = genome[chrom] - lengths

        # same scaffold + same length as every real deletion
        rstart = np.floor(rng.random(len(lengths)) * (maxstart + 1)).astype(np.int64)
        rend = rstart + lengths

        for i, cls in enumerate(classes):
            ov = overlap_bp(rstart, rend, models[cls].get(chrom))
            null_bp[perm,i] += ov.sum()
            null_events[perm,i] += np.sum(ov > 0)

    if (perm + 1) % 100 == 0:
        print(f"Completed {perm+1}/{a.nperm} permutations")

# ---------- Benjamini-Hochberg ----------
def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * len(p) / np.arange(1, len(p)+1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    result = np.empty_like(q)
    result[order] = np.minimum(q, 1)
    return result

p_bp = np.array([(1 + np.sum(null_bp[:,i] >= obs_bp[i])) / (a.nperm + 1) for i in range(len(classes))])
p_ev = np.array([(1 + np.sum(null_events[:,i] >= obs_events[i])) / (a.nperm + 1) for i in range(len(classes))])
q_bp, q_ev = bh(p_bp), bh(p_ev)

# ---------- summary ----------
outfile = os.path.join(a.outdir, "deletion_repeat_enrichment.tsv")

with open(outfile, "w") as o:
    o.write("CLASS\tOBS_EVENTS\tOBS_EVENT_PCT\tNULL_EVENT_PCT\tEVENT_FOLD\tP_EVENT\tFDR_EVENT\tOBS_BP\tOBS_BP_PCT\tNULL_BP_PCT\tBP_FOLD\tP_BP\tFDR_BP\n")

    for i, cls in enumerate(classes):
        obs_ep = 100 * obs_events[i] / total_events
        null_ep = 100 * np.mean(null_events[:,i]) / total_events
        obs_bpp = 100 * obs_bp[i] / total_bp
        null_bpp = 100 * np.mean(null_bp[:,i]) / total_bp

        event_fold = obs_ep / null_ep if null_ep else np.nan
        bp_fold = obs_bpp / null_bpp if null_bpp else np.nan

        o.write(
            f"{cls}\t{obs_events[i]}\t{obs_ep:.3f}\t{null_ep:.3f}\t{event_fold:.3f}\t{p_ev[i]:.6g}\t{q_ev[i]:.6g}\t"
            f"{obs_bp[i]}\t{obs_bpp:.3f}\t{null_bpp:.3f}\t{bp_fold:.3f}\t{p_bp[i]:.6g}\t{q_bp[i]:.6g}\n"
        )

# Save permutation distributions for reproducibility.
nullfile = os.path.join(a.outdir, "deletion_repeat_null.tsv")
with open(nullfile, "w") as o:
    header = ["PERM"] + [f"{c}_BP_PCT" for c in classes]
    o.write("\t".join(header) + "\n")
    for n in range(a.nperm):
        pct = 100 * null_bp[n] / total_bp
        o.write(str(n+1) + "\t" + "\t".join(f"{x:.5f}" for x in pct) + "\n")

print("\nDone.")
print("Summary:", outfile)
print("Null distributions:", nullfile)