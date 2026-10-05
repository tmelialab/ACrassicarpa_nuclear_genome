#!/bin/bash
#SBATCH --job-name=acra_genecontent
#SBATCH --output=acra_genecontent_%j.out
#SBATCH --error=acra_genecontent_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

# Usage:
# sbatch gene_content_orthofinder.slurm \
#   V4.gtf V4.proteins.fa \
#   Acra3RX.gtf Acra3RX.proteins.fa \
#   OUTDIR

V4_GTF="$1"
V4_PROT="$2"
REF_GTF="$3"
REF_PROT="$4"
OUTDIR="$5"

mkdir -p "$OUTDIR/proteomes"

# ============================================================
# Activate OrthoFinder environment
# ============================================================

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate orthofinder

echo "OrthoFinder:"
which orthofinder
orthofinder --version || true

V4_PRIMARY="$OUTDIR/proteomes/v4.fa"
REF_PRIMARY="$OUTDIR/proteomes/Acra3RX.fa"

# ============================================================
# Select longest protein isoform per gene
# ============================================================

python3 - \
"$V4_GTF" "$V4_PROT" "$V4_PRIMARY" \
"$REF_GTF" "$REF_PROT" "$REF_PRIMARY" <<'PY'

import sys, re

v4_gtf, v4_fa, v4_out, ref_gtf, ref_fa, ref_out = sys.argv[1:]

def read_fasta(fn):
    seqs = {}
    name = None
    buf = []

    with open(fn) as f:
        for line in f:
            line = line.rstrip()

            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)

                name = line[1:].split()[0]
                buf = []

            else:
                buf.append(line)

        if name is not None:
            seqs[name] = "".join(buf)

    return seqs


def transcript_gene_map(gtf):
    m = {}

    with open(gtf) as f:
        for line in f:

            if not line.strip() or line.startswith("#"):
                continue

            x = line.rstrip().split(None,8)

            if len(x) < 9:
                continue

            feature = x[2]
            attr = x[8]

            tid = None
            gid = None

            z = re.search(r'transcript_id\s+"([^"]+)"', attr)
            if z:
                tid = z.group(1)

            z = re.search(r'gene_id\s+"([^"]+)"', attr)
            if z:
                gid = z.group(1)

            # BRAKER transcript lines may contain only g1.t1
            if feature in ("transcript","mRNA") and tid is None:
                token = attr.strip().strip(";")
                if token and " " not in token:
                    tid = token

            if tid and gid:
                m[tid] = gid

            elif tid and ".t" in tid:
                m[tid] = tid.rsplit(".t",1)[0]

    return m


def write_longest(gtf, fasta, outfile, label):

    seqs = read_fasta(fasta)
    t2g = transcript_gene_map(gtf)

    best = {}
    unmapped = 0

    for tid,seq in seqs.items():

        gid = t2g.get(tid)

        if gid is None and ".t" in tid:
            gid = tid.rsplit(".t",1)[0]

        if gid is None:
            unmapped += 1
            continue

        if gid not in best or len(seq) > len(best[gid][1]):
            best[gid] = (tid,seq)

    with open(outfile,"w") as o:
        for gid,(tid,seq) in sorted(best.items()):
            o.write(f">{gid}\n")

            for i in range(0,len(seq),60):
                o.write(seq[i:i+60]+"\n")

    print(label)
    print("  input proteins:", len(seqs))
    print("  genes represented:", len(best))
    print("  unmapped:", unmapped)


write_longest(v4_gtf,v4_fa,v4_out,"v4")
write_longest(ref_gtf,ref_fa,ref_out,"Acra3RX")

PY

echo
echo "======================================"
echo "Representative protein sets"
echo "======================================"

echo -n "v4 genes:      "
grep -c '^>' "$V4_PRIMARY"

echo -n "Acra3RX genes: "
grep -c '^>' "$REF_PRIMARY"

# ============================================================
# OrthoFinder
# ============================================================

orthofinder \
    -f "$OUTDIR/proteomes" \
    -t ${SLURM_CPUS_PER_TASK} \
    -a ${SLURM_CPUS_PER_TASK}

echo
echo "======================================"
echo "OrthoFinder completed"
echo "======================================"

RESULT=$(find "$OUTDIR/proteomes" \
    -type d -name 'Results_*' | sort | tail -1)

echo "Result directory:"
echo "$RESULT"

echo
echo "Overall statistics:"
cat "$RESULT/Comparative_Genomics_Statistics/Statistics_Overall.tsv" 2>/dev/null || true