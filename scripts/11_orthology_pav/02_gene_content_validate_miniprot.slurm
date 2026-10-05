#!/bin/bash
#SBATCH --job-name=genecontent_val
#SBATCH --output=genecontent_val_%j.out
#SBATCH --error=genecontent_val_%j.err
#SBATCH --cpus-per-task=8

set -euo pipefail

# Validate candidate gene-content differences using miniprot
# and associate residual candidates with supported SVs.
#
# Usage:
# sbatch 02_gene_content_validate_miniprot.slurm \
#   <orthofinder_results> \
#   <AcraUNRI.fa> \
#   <Acra3RX_reference.fa> \
#   <Acra3RX_haplotype1.fa> \
#   <Acra3RX_haplotype2.fa> \
#   <output_dir> \
#   <Acra3RX.gtf> \
#   <supported_SVs.vcf> \
#   <SVIM_assembly_SVs.vcf>

RESULTS="${1:?Missing OrthoFinder results directory}"
V4="${2:?Missing AcraUNRI genome}"
REF="${3:?Missing Acra3RX reference genome}"
HAP1="${4:?Missing Acra3RX haplotype 1}"
HAP2="${5:?Missing Acra3RX haplotype 2}"
OUT="${6:?Missing output directory}"
REF_GTF="${7:?Missing Acra3RX GTF}"
SUPPORTED_VCF="${8:?Missing supported SV VCF}"
SVIM_VCF="${9:?Missing SVIM-asm source VCF}"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
MINIPROT="${MINIPROT:-miniprot}"

command -v "$MINIPROT" >/dev/null 2>&1 || {
    echo "ERROR: miniprot not found." >&2
    exit 1
}

mkdir -p "$OUT"

PROTDIR=$(dirname "$(dirname "$RESULTS")")
V4_PROT="$PROTDIR/v4.fa"
REF_PROT="$PROTDIR/Acra3RX.fa"

OG="$RESULTS/Orthogroups/Orthogroups.tsv"
UN="$RESULTS/Orthogroups/Orthogroups_UnassignedGenes.tsv"


# ============================================================
# 1. Extract OrthoFinder candidates
# ============================================================

python3 - "$OG" "$UN" "$V4_PROT" "$REF_PROT" "$OUT" <<'PY'
import sys,csv
og,un,v4fa,reffa,out=sys.argv[1:]

def fasta(fn):
    d={}; n=None
    for l in open(fn):
        l=l.strip()
        if l.startswith(">"):
            n=l[1:].split()[0]; d[n]=""
        elif n: d[n]+=l
    return d

def genes(x):
    return x.replace(","," ").split() if x else []

v4=set(); ref=set()

with open(og) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        a=genes(r.get("v4","")); b=genes(r.get("Acra3RX",""))
        if a and not b: v4.update(a)
        if b and not a: ref.update(b)

with open(un) as f:
    for r in csv.DictReader(f,delimiter="\t"):
        v4.update(genes(r.get("v4","")))
        ref.update(genes(r.get("Acra3RX","")))

def write(ids,seqs,fn):
    with open(fn,"w") as o:
        for g in sorted(ids):
            if g not in seqs: continue
            s=seqs[g]
            o.write(f">{g}\n")
            for i in range(0,len(s),60):
                o.write(s[i:i+60]+"\n")

write(v4,fasta(v4fa),out+"/v4_candidates.fa")
write(ref,fasta(reffa),out+"/Acra3RX_candidates.fa")

print("v4 candidates:",len(v4))
print("Acra3RX candidates:",len(ref))
PY


# ============================================================
# 2. Miniprot indexes + mappings
# ============================================================

index () {
    [ -s "$2" ] || "$MINIPROT" -t "$THREADS" -d "$2" "$1"
}

mapprot () {
    [ -s "$3" ] || \
    "$MINIPROT" -I -t "$THREADS" --gff-only --outn 5 \
        "$1" "$2" > "$3"
}

index "$V4"   "$OUT/v4.mpi"
index "$REF"  "$OUT/ref.mpi"
index "$HAP1" "$OUT/hap1.mpi"
index "$HAP2" "$OUT/hap2.mpi"

mapprot "$OUT/ref.mpi"  "$OUT/v4_candidates.fa" "$OUT/v4_to_ref.gff3"
mapprot "$OUT/hap1.mpi" "$OUT/v4_candidates.fa" "$OUT/v4_to_hap1.gff3"
mapprot "$OUT/hap2.mpi" "$OUT/v4_candidates.fa" "$OUT/v4_to_hap2.gff3"
mapprot "$OUT/v4.mpi" "$OUT/Acra3RX_candidates.fa" "$OUT/ref_to_v4.gff3"


# ============================================================
# 3. Final gene-content classification
# ============================================================

python3 - \
"$OUT/v4_candidates.fa" "$OUT/Acra3RX_candidates.fa" \
"$OUT/v4_to_ref.gff3" "$OUT/v4_to_hap1.gff3" \
"$OUT/v4_to_hap2.gff3" "$OUT/ref_to_v4.gff3" "$OUT" <<'PY'

import sys,re
from collections import Counter

v4fa,reffa,g0,g1,g2,gr,out=sys.argv[1:]

def fasta(fn):
    d={}; n=None
    for l in open(fn):
        l=l.strip()
        if l.startswith(">"):
            n=l[1:].split()[0]; d[n]=""
        elif n: d[n]+=l
    return d

def attrs(s):
    return dict(x.split("=",1) for x in s.split(";") if "=" in x)

def tier(cov,ident,fs,stop):
    if cov>=.80 and ident>=.90 and fs==0 and stop==0: return 2
    if cov>=.30 and ident>=.50: return 1
    return 0

def parse(fn,seqs,source):
    H={}
    for l in open(fn):
        if l.startswith("#"): continue
        x=l.rstrip().split("\t")
        if len(x)<9 or x[2]!="mRNA": continue

        a=attrs(x[8])
        m=re.match(r'(\S+)\s+(\d+)\s+(\d+)',a.get("Target",""))
        if not m or m.group(1) not in seqs: continue

        g=m.group(1)
        cov=(abs(int(m.group(3))-int(m.group(2)))+1)/len(seqs[g])
        ident=float(a.get("Identity",0))
        if ident>1: ident/=100

        fs=int(a.get("Frameshift",0))
        stop=int(a.get("StopCodon",0))

        h=(tier(cov,ident,fs,stop),cov,ident,-fs,-stop,source,fs,stop)
        if g not in H or h[:5]>H[g][:5]:
            H[g]=h
    return H

def combine(*sets):
    R={}
    for D in sets:
        for g,h in D.items():
            if g not in R or h[:5]>R[g][:5]:
                R[g]=h
    return R

def report(name,seqs,H):
    c=Counter(); final=[]

    with open(out+"/"+name+"_validation.tsv","w") as o:
        o.write("gene\tclass\tsource\tcoverage\tidentity\tframeshift\tstop\n")

        for g in sorted(seqs):
            h=H.get(g)

            if not h or h[0]==0:
                cls="no_convincing_counterpart"
            elif h[0]==2:
                cls="intact_like_counterpart"
            else:
                cls="partial_or_divergent_counterpart"

            c[cls]+=1

            if cls=="no_convincing_counterpart":
                final.append(g)

            if h:
                o.write(f"{g}\t{cls}\t{h[5]}\t{h[1]:.4f}\t"
                        f"{h[2]:.4f}\t{h[6]}\t{h[7]}\n")
            else:
                o.write(f"{g}\t{cls}\t.\t0\t0\t0\t0\n")

    with open(out+"/"+name+"_final_PAV.txt","w") as o:
        for g in final: o.write(g+"\n")

    with open(out+"/"+name+"_final_PAV.fa","w") as o:
        for g in final:
            o.write(f">{g}\n")
            s=seqs[g]
            for i in range(0,len(s),60):
                o.write(s[i:i+60]+"\n")

    return c,len(final)

v4=fasta(v4fa); ref=fasta(reffa)

Hv4=combine(
    parse(g0,v4,"Acra3RX_100k"),
    parse(g1,v4,"Acra3RX_hap1"),
    parse(g2,v4,"Acra3RX_hap2")
)

Href=parse(gr,ref,"v4")

c1,n1=report("v4",v4,Hv4)
c2,n2=report("Acra3RX",ref,Href)

with open(out+"/gene_content_final_summary.tsv","w") as o:
    o.write("genome\ttotal\tintact_like\tpartial_divergent\tno_counterpart\n")
    o.write(f"v4\t{len(v4)}\t{c1['intact_like_counterpart']}\t"
            f"{c1['partial_or_divergent_counterpart']}\t{n1}\n")
    o.write(f"Acra3RX\t{len(ref)}\t{c2['intact_like_counterpart']}\t"
            f"{c2['partial_or_divergent_counterpart']}\t{n2}\n")
PY


# ============================================================
# 4. Acra3RX candidate genes vs supported DELs
# ============================================================

python3 - "$REF_GTF" "$SUPPORTED_VCF" \
"$OUT/Acra3RX_final_PAV.txt" "$OUT/Acra3RX_PAV_DEL.tsv" <<'PY'

import sys,re,gzip
gtf,vcf,ids,out=sys.argv[1:]

wanted=set(x.strip() for x in open(ids) if x.strip())
genes={}

for l in open(gtf):
    if not l.strip() or l.startswith("#"): continue
    x=l.rstrip().split(None,8)
    if len(x)<9 or x[2]!="gene": continue

    a=x[8]
    m=re.search(r'gene_id\s+"([^"]+)"',a)
    if m: g=m.group(1)
    else:
        m=re.search(r'ID=([^;]+)',a)
        g=m.group(1) if m else a.strip().split(";")[0].strip('"')

    if g in wanted:
        genes[g]=(x[0],int(x[3]),int(x[4]))

op=gzip.open if vcf.endswith(".gz") else open
dels={}

with op(vcf,"rt") as f:
    for l in f:
        if l.startswith("#"): continue
        x=l.rstrip().split("\t")
        D=dict(z.split("=",1) for z in x[7].split(";") if "=" in z)

        if D.get("SVTYPE")!="DEL": continue

        s=int(x[1])
        L=abs(int(D.get("SVLEN","0").split(",")[0]))
        e=int(D.get("END",s+L-1))

        dels.setdefault(x[0],[]).append((s,e,x[2]))

with open(out,"w") as o:
    o.write("gene\tchrom\tgene_start\tgene_end\t"
            "del_start\tdel_end\tdel_id\trelation\n")

    for g,(chrom,s,e) in genes.items():
        for ds,de,did in dels.get(chrom,[]):
            if ds<=e and de>=s:
                rel="whole_gene" if ds<=s and de>=e else "partial_overlap"
                o.write(f"{g}\t{chrom}\t{s}\t{e}\t"
                        f"{ds}\t{de}\t{did}\t{rel}\n")
PY


# ============================================================
# 5. Recover sequence for each supported INS
# ============================================================

python3 - "$SUPPORTED_VCF" "$SVIM_VCF" \
"$OUT/supported_insertions.fa" "$OUT/supported_insertions.tsv" <<'PY'

import sys,gzip

supported_file,source_file,fa,meta=sys.argv[1:]

def op(fn):
    return gzip.open(fn,"rt") if fn.endswith(".gz") else open(fn)

def info(s):
    return dict(x.split("=",1) for x in s.split(";") if "=" in x)

def getseq(ref,alt,I):
    for k in ("SEQ","SEQS","INSSEQ","SVINSSEQ","CONSENSUS"):
        v=I.get(k)
        if v and v!=".":
            s=v.split(",")[0]
            if not s.startswith("<"):
                return s

    if alt and alt!="." and not alt.startswith("<") and len(alt)>len(ref):
        return alt[len(ref):] if alt.startswith(ref) else alt

    return None

def read_ins(fn):
    R=[]
    with op(fn) as f:
        for l in f:
            if l.startswith("#"): continue

            x=l.rstrip().split("\t")
            I=info(x[7])

            if I.get("SVTYPE")!="INS":
                continue

            try:
                L=abs(int(I.get("SVLEN","0").split(",")[0]))
            except:
                L=0

            R.append({
                "chrom":x[0],
                "pos":int(x[1]),
                "id":x[2],
                "len":L,
                "ref":x[3],
                "alt":x[4],
                "info":I,
                "seq":getseq(x[3],x[4],I)
            })
    return R

supported=read_ins(supported_file)
source=read_ins(source_file)

by_id={}
by_key={}

for r in source:
    if r["id"] not in ("","."):
        by_id[r["id"]]=r
    by_key.setdefault((r["chrom"],r["pos"],r["len"]),[]).append(r)

resolved=[]

for r in supported:

    seq=r["seq"]
    src_id=r["id"]

    # First try exact VCF ID
    if not seq and r["id"] not in ("",".") and r["id"] in by_id:
        q=by_id[r["id"]]
        seq=q["seq"]
        src_id=q["id"]

    # Fallback: chromosome + position + SVLEN
    if not seq:
        candidates=by_key.get((r["chrom"],r["pos"],r["len"]),[])
        candidates=[q for q in candidates if q["seq"]]

        if candidates:
            q=min(candidates,key=lambda z:abs(len(z["seq"])-r["len"]))
            seq=q["seq"]
            src_id=q["id"]

    if seq:
        resolved.append((r,src_id,seq))

with open(fa,"w") as F, open(meta,"w") as M:

    M.write("ins_id\tchrom\tpos\tsvlen\tsupported_vcf_id\tsource_vcf_id\n")

    for i,(r,src_id,seq) in enumerate(resolved,1):

        iid=f"SUPP_INS{i:06d}"

        F.write(f">{iid}\n")
        for j in range(0,len(seq),80):
            F.write(seq[j:j+80]+"\n")

        M.write(
            f"{iid}\t{r['chrom']}\t{r['pos']}\t{r['len']}\t"
            f"{r['id']}\t{src_id}\n"
        )

print("ONT-supported INS:",len(supported))
print("Sequence-resolved INS recovered:",len(resolved))
PY


# ============================================================
# 6. v4 candidate PAV proteins vs supported INS sequence
# ============================================================

rm -f "$OUT/supported_INS.mpi"
rm -f "$OUT/v4_PAV_to_supported_INS.gff3"
rm -f "$OUT/v4_PAV_INS.tsv"

if [ -s "$OUT/supported_insertions.fa" ] && \
   [ -s "$OUT/v4_final_PAV.fa" ]; then

    "$MINIPROT" -t "$THREADS" \
        -d "$OUT/supported_INS.mpi" \
        "$OUT/supported_insertions.fa"

    "$MINIPROT" -I -t "$THREADS" --gff-only --outn 5 \
        "$OUT/supported_INS.mpi" \
        "$OUT/v4_final_PAV.fa" \
        > "$OUT/v4_PAV_to_supported_INS.gff3"

    python3 - \
    "$OUT/v4_final_PAV.fa" \
    "$OUT/v4_PAV_to_supported_INS.gff3" \
    "$OUT/supported_insertions.tsv" \
    "$OUT/v4_PAV_INS.tsv" <<'PY'

import sys,re

fa,gff,meta,out=sys.argv[1:]

seqs={}; n=None

for l in open(fa):
    l=l.strip()
    if l.startswith(">"):
        n=l[1:].split()[0]; seqs[n]=""
    elif n:
        seqs[n]+=l

M={}

for l in open(meta):
    if l.startswith("ins_id"): continue
    x=l.rstrip().split("\t")
    M[x[0]]=x[1:]

best={}

for l in open(gff):

    if l.startswith("#"): continue

    x=l.rstrip().split("\t")

    if len(x)<9 or x[2]!="mRNA":
        continue

    A=dict(z.split("=",1) for z in x[8].split(";") if "=" in z)

    m=re.match(r'(\S+)\s+(\d+)\s+(\d+)',A.get("Target",""))

    if not m or m.group(1) not in seqs:
        continue

    g=m.group(1)

    cov=(abs(int(m.group(3))-int(m.group(2)))+1)/len(seqs[g])

    ident=float(A.get("Identity",0))
    if ident>1: ident/=100

    fs=int(A.get("Frameshift",0))
    stop=int(A.get("StopCodon",0))

    tier=2 if cov>=.80 and ident>=.90 and fs==0 and stop==0 else \
         1 if cov>=.30 and ident>=.50 else 0

    h=(tier,cov,ident,-fs,-stop,x[0])

    if g not in best or h[:5]>best[g][:5]:
        best[g]=h

with open(out,"w") as o:

    o.write("gene\tclass\tins_id\tchrom\tpos\tsvlen\tcoverage\tidentity\n")

    for g in sorted(seqs):

        h=best.get(g)

        if not h or h[0]==0:
            o.write(f"{g}\tno_INS_association\t.\t.\t.\t.\t0\t0\n")
            continue

        cls="intact_like_INS" if h[0]==2 else "partial_INS"

        iid=h[5]
        chrom,pos,L,support_id,source_id=M.get(
            iid,[".",".",".",".","."]
        )

        o.write(
            f"{g}\t{cls}\t{iid}\t{chrom}\t{pos}\t{L}\t"
            f"{h[1]:.4f}\t{h[2]:.4f}\n"
        )
PY

fi


# ============================================================
# 7. Final summary
# ============================================================

python3 - "$OUT" <<'PY'

import sys,os
out=sys.argv[1:]

gc={}

for l in open(out+"/gene_content_final_summary.tsv"):
    if l.startswith("genome"): continue
    x=l.rstrip().split("\t")
    gc[x[0]]=x

delgenes=set()
whole=set()

fn=out+"/Acra3RX_PAV_DEL.tsv"

if os.path.exists(fn):
    for l in open(fn):
        if l.startswith("gene"): continue
        x=l.rstrip().split("\t")
        delgenes.add(x[0])

        if x[7]=="whole_gene":
            whole.add(x[0])

ins_full=set()
ins_partial=set()

fn=out+"/v4_PAV_INS.tsv"

if os.path.exists(fn):
    for l in open(fn):
        if l.startswith("gene"): continue

        x=l.rstrip().split("\t")

        if x[1]=="intact_like_INS":
            ins_full.add(x[0])

        elif x[1]=="partial_INS":
            ins_partial.add(x[0])

with open(out+"/gene_content_SV_summary.tsv","w") as o:

    o.write("metric\tcount\n")

    o.write(f"v4_final_PAV\t{gc['v4'][4]}\n")
    o.write(f"v4_PAV_intact_INS_association\t{len(ins_full)}\n")
    o.write(f"v4_PAV_partial_INS_association\t{len(ins_partial)}\n")

    o.write(f"Acra3RX_final_PAV\t{gc['Acra3RX'][4]}\n")
    o.write(f"Acra3RX_PAV_overlapping_DEL\t{len(delgenes)}\n")
    o.write(f"Acra3RX_PAV_whole_gene_DEL\t{len(whole)}\n")

print("\nGene-content summary:")
print(open(out+"/gene_content_final_summary.tsv").read())

print("SV-association summary:")
print(open(out+"/gene_content_SV_summary.tsv").read())

PY

echo "Done."