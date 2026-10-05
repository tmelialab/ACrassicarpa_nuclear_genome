#!/usr/bin/env python3

import argparse, gzip, os
from collections import defaultdict, Counter
from bisect import bisect_left

p = argparse.ArgumentParser()
p.add_argument("--vcf", required=True, help="ONT-supported VCF")
p.add_argument("--gtf", required=True, help="Acra3RX GTF")
p.add_argument("--outdir", default="output")
p.add_argument("--distance", type=int, default=5000, help="near-gene distance (bp)")
a = p.parse_args()

os.makedirs(a.outdir, exist_ok=True)

def opener(x):
    return gzip.open(x, "rt") if x.endswith(".gz") else open(x)

def gene_id(attr):
    import re
    m = re.search(r'gene_id\s+"([^"]+)"', attr)
    if m: return m.group(1)
    m = re.search(r'gene_id[=\s]+([^; "\t]+)', attr)
    if m: return m.group(1)
    return attr.strip().split()[0].strip('";') if attr.strip() else "NA"

# ============================================================
# 1. Read GTF
# ============================================================

genes = defaultdict(list)
exons = defaultdict(list)
cds = defaultdict(list)

with open(a.gtf) as f:
    for line in f:
        if not line.strip() or line.startswith("#"): continue

        # Works for both tab- and whitespace-delimited GTF
        x = line.rstrip().split(None, 8)
        if len(x) < 9: continue

        chrom, feature = x[0], x[2]
        start, end = int(x[3])-1, int(x[4])       # 0-based half-open
        gid = gene_id(x[8])

        if feature == "gene": genes[chrom].append((start,end,gid))
        elif feature == "exon": exons[chrom].append((start,end,gid))
        elif feature == "CDS": cds[chrom].append((start,end,gid))

for D in (genes, exons, cds):
    for chrom in D: D[chrom].sort()

print("Genes:", sum(len(x) for x in genes.values()))
print("Exons:", sum(len(x) for x in exons.values()))
print("CDS:", sum(len(x) for x in cds.values()))

# ============================================================
# 2. Interval helpers
# ============================================================

def overlaps(intervals, start, end):
    """Return unique gene IDs overlapping [start,end)."""
    if not intervals: return set()
    starts = [x[0] for x in intervals]
    i = bisect_left(starts, end)
    out = set()
    for s,e,g in intervals[:i]:
        if e > start: out.add(g)
    return out

def full_genes(intervals, start, end):
    return {g for s,e,g in intervals if s >= start and e <= end}

# Precompute information for nearest-gene lookup
nearest_index = {}

for chrom, g in genes.items():
    starts = [x[0] for x in g]
    prefix_end, prefix_gene = [], []
    best_end, best_gene = -1, None

    for s,e,gid in g:
        if e > best_end: best_end, best_gene = e, gid
        prefix_end.append(best_end); prefix_gene.append(best_gene)

    nearest_index[chrom] = (g, starts, prefix_end, prefix_gene)

def nearest_gene(chrom, start, end):
    if chrom not in nearest_index: return None, None
    g, starts, pend, pgene = nearest_index[chrom]
    i = bisect_left(starts, end)

    candidates = []

    if i < len(g):
        dist = max(0, g[i][0] - end)
        candidates.append((dist, g[i][2]))

    if i > 0:
        dist = max(0, start - pend[i-1])
        candidates.append((dist, pgene[i-1]))

    if not candidates: return None, None
    return min(candidates)

# ============================================================
# 3. Read VCF and classify SVs
# ============================================================

events = []
relations = []
gene_events = defaultdict(lambda: {"DEL":set(), "INS":set(), "impact":set()})

def parse_info(s):
    d = {}
    for z in s.split(";"):
        if "=" in z:
            k,v = z.split("=",1); d[k] = v
    return d

severity = {
    "whole_gene_deletion": 1,
    "CDS_overlap": 2,
    "exon_overlap": 3,
    "intronic_overlap": 4,
    "CDS_insertion": 2,
    "exon_insertion": 3,
    "intronic_insertion": 4,
    "near_gene": 5,
    "intergenic": 6
}

with opener(a.vcf) as f:
    for n,line in enumerate(f,1):
        if line.startswith("#"): continue
        x = line.rstrip().split("\t")
        chrom, pos, vid = x[0], int(x[1]), x[2]
        info = parse_info(x[7])
        typ = info.get("SVTYPE","")
        if typ not in ("DEL","INS"): continue

        svlen = abs(int(info.get("SVLEN","0").split(",")[0]))
        sid = vid if vid != "." else f"{typ}_{chrom}_{pos}_{n}"

        fmt = x[8].split(":"); sample = x[9].split(":")
        gt = sample[fmt.index("GT")] if "GT" in fmt and fmt.index("GT") < len(sample) else "NA"

        impacts, affected = {}, set()

        # ----------------------------------------------------
        # DEL: actual deleted interval excludes VCF anchor base
        # ----------------------------------------------------
        if typ == "DEL":
            end = int(info.get("END", pos + svlen))
            start = pos

            gset = overlaps(genes.get(chrom,[]), start, end)
            eset = overlaps(exons.get(chrom,[]), start, end)
            cset = overlaps(cds.get(chrom,[]), start, end)
            fset = full_genes(genes.get(chrom,[]), start, end)

            for gid in gset:
                if gid in fset: impact = "whole_gene_deletion"
                elif gid in cset: impact = "CDS_overlap"
                elif gid in eset: impact = "exon_overlap"
                else: impact = "intronic_overlap"

                impacts[gid] = impact
                affected.add(gid)
                relations.append((sid,typ,chrom,pos,svlen,gid,impact))

        # ----------------------------------------------------
        # INS: breakpoint immediately after VCF POS
        # ----------------------------------------------------
        else:
            start = pos
            end = pos + 1

            gset = overlaps(genes.get(chrom,[]), start, end)
            eset = overlaps(exons.get(chrom,[]), start, end)
            cset = overlaps(cds.get(chrom,[]), start, end)

            for gid in gset:
                if gid in cset: impact = "CDS_insertion"
                elif gid in eset: impact = "exon_insertion"
                else: impact = "intronic_insertion"

                impacts[gid] = impact
                affected.add(gid)
                relations.append((sid,typ,chrom,pos,svlen,gid,impact))

        # ----------------------------------------------------
        # No direct gene overlap: classify proximity
        # ----------------------------------------------------
        if not affected:
            dist, gid = nearest_gene(chrom,start,end)

            if dist is not None and dist <= a.distance:
                main = "near_gene"
                gene_list = gid
                relations.append((sid,typ,chrom,pos,svlen,gid,f"near_gene_{dist}bp"))
                gene_events[gid][typ].add(sid)
                gene_events[gid]["impact"].add("near_gene")
            else:
                main = "intergenic"
                gene_list = "."
        else:
            main = min(impacts.values(), key=lambda z: severity[z])
            gene_list = ",".join(sorted(affected))

            for gid,impact in impacts.items():
                gene_events[gid][typ].add(sid)
                gene_events[gid]["impact"].add(impact)

        events.append((sid,chrom,pos,typ,svlen,gt,main,len(affected),gene_list))

# ============================================================
# 4. Event-level output
# ============================================================

event_file = os.path.join(a.outdir,"sv_gene_impact.tsv")
with open(event_file,"w") as o:
    o.write("SV_ID\tCHROM\tPOS\tSVTYPE\tSVLEN\tGT\tIMPACT\tN_OVERLAPPED_GENES\tGENE_IDS\n")
    for r in events: o.write("\t".join(map(str,r))+"\n")

# ============================================================
# 5. SV-gene relations
# ============================================================

relation_file = os.path.join(a.outdir,"sv_gene_relations.tsv")
with open(relation_file,"w") as o:
    o.write("SV_ID\tSVTYPE\tCHROM\tPOS\tSVLEN\tGENE_ID\tIMPACT\n")
    for r in relations: o.write("\t".join(map(str,r))+"\n")

# ============================================================
# 6. Summary by SV type and impact
# ============================================================

counts = Counter((r[3],r[6]) for r in events)
totals = Counter(r[3] for r in events)

summary_file = os.path.join(a.outdir,"sv_gene_impact_summary.tsv")
with open(summary_file,"w") as o:
    o.write("SVTYPE\tIMPACT\tCOUNT\tPERCENT\n")
    for typ in ("DEL","INS"):
        for impact in sorted({k[1] for k in counts if k[0] == typ}, key=lambda x: severity[x]):
            n = counts[(typ,impact)]
            o.write(f"{typ}\t{impact}\t{n}\t{100*n/totals[typ]:.2f}\n")

# ============================================================
# 7. Unique genes associated with SVs
# ============================================================

gene_file = os.path.join(a.outdir,"sv_affected_genes.tsv")
with open(gene_file,"w") as o:
    o.write("GENE_ID\tDEL_EVENTS\tINS_EVENTS\tTOTAL_EVENTS\tIMPACTS\n")
    for gid,d in sorted(gene_events.items()):
        nd, ni = len(d["DEL"]), len(d["INS"])
        o.write(f"{gid}\t{nd}\t{ni}\t{nd+ni}\t{','.join(sorted(d['impact']))}\n")

# ============================================================
# 8. Overall overview
# ============================================================

direct = sum(r[6] not in ("near_gene","intergenic") for r in events)
near = sum(r[6] == "near_gene" for r in events)
inter = sum(r[6] == "intergenic" for r in events)

overview_file = os.path.join(a.outdir,"sv_gene_overview.txt")
with open(overview_file,"w") as o:
    o.write(f"Total INS/DEL SVs: {len(events)}\n")
    o.write(f"DEL: {totals['DEL']}\n")
    o.write(f"INS: {totals['INS']}\n")
    o.write(f"Direct gene overlap: {direct} ({100*direct/len(events):.2f}%)\n")
    o.write(f"Within {a.distance} bp of gene: {near} ({100*near/len(events):.2f}%)\n")
    o.write(f"Intergenic: {inter} ({100*inter/len(events):.2f}%)\n")
    o.write(f"Unique genes associated with SVs: {len(gene_events)}\n")

print("\n==============================")
print("SV–gene impact analysis")
print("==============================")
print(open(overview_file).read())
print("Impact summary:")
print(open(summary_file).read())