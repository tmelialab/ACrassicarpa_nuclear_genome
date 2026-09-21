#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

def parse_gtf_attrs(attr):
    """
    Parse standard GTF attributes: key "value";
    Returns dict.
    """
    d = {}
    # matches: key "value"
    for m in re.finditer(r'(\S+)\s+"([^"]+)"', attr):
        d[m.group(1)] = m.group(2)
    return d

def read_gene_transcript_map(gtf_path):
    """
    Returns:
      tx2gene: transcript_id -> gene_id
    Supports AUGUSTUS/BRAKER GTF variants where:
      - gene/transcript IDs may be in column 9 as plain IDs (e.g. 'g1', 'g1.t1')
      - or present as transcript_id/gene_id attributes
    """
    tx2gene = {}

    gene_id_current = None
    transcript_id_current = None

    with open(gtf_path, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            attr = cols[8].strip()

            attrs = parse_gtf_attrs(attr)

            if feature == "gene":
                # AUGUSTUS sometimes stores gene id as plain 'g1' in column 9
                # and may not include gene_id attr
                gene_id_current = attrs.get("gene_id", attr.split()[0] if attr else None)

            elif feature in ("transcript", "mRNA"):
                transcript_id_current = attrs.get("transcript_id", attr.split()[0] if attr else None)
                # gene id may be in attrs or from last seen gene feature
                gid = attrs.get("gene_id", gene_id_current)
                tid = transcript_id_current
                if gid and tid:
                    tx2gene[tid] = gid

            else:
                # For CDS/exon/stop_codon lines, they usually contain transcript_id and gene_id
                tid = attrs.get("transcript_id", None)
                gid = attrs.get("gene_id", None)
                if tid and gid:
                    tx2gene[tid] = gid

    return tx2gene

def read_fasta_ids(fa_path):
    ids = []
    with open(fa_path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids

def read_diamond_besthit(tsv_path):
    """
    Reads DIAMOND outfmt 6 TSV where qseqid is col1.
    Keeps best bitscore per query.
    Returns: dict q -> (sseqid, evalue, bitscore, stitle)
    """
    hits = {}
    with open(tsv_path) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            q = cols[0]
            s = cols[1]
            if len(cols) >= 12:
                evalue = cols[10]
                bitscore = cols[11]
                stitle = cols[12] if len(cols) > 12 else s
            else:
                evalue = cols[4]
                bitscore = cols[5]
                stitle = cols[6] if len(cols) > 6 else s
            if q not in hits:
                hits[q] = (s, evalue, bitscore, stitle)
            else:
                try:
                    if float(bitscore) > float(hits[q][2]):
                        hits[q] = (s, evalue, bitscore, stitle)
                except:
                    pass
    return hits

def read_interproscan_tsv(ipr_path):
    """
    Aggregates Pfam IDs, InterPro IDs, GO terms, and descriptions per protein.
    Assumes InterProScan TSV with >=14 columns, with:
      1: protein_id
      4: analysis
      5: signature accession
      6: signature description
      12: interpro accession
      13: interpro description
      14: GO terms
    (This matches typical InterProScan TSV; if your column positions differ, we can adjust.)
    """
    pfam_ids = defaultdict(set)
    ipr_ids  = defaultdict(set)
    go_terms = defaultdict(set)
    sig_desc = defaultdict(set)
    ipr_desc = defaultdict(set)

    with open(ipr_path) as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 14:
                continue
            pid = cols[0]
            analysis = cols[3]
            sig_acc  = cols[4]
            sig_d    = cols[5]
            ipr_acc  = cols[11]
            ipr_d    = cols[12]
            go       = cols[13]

            if analysis.lower() == "pfam" and sig_acc and sig_acc != "-":
                pfam_ids[pid].add(sig_acc)
            if sig_d and sig_d != "-":
                sig_desc[pid].add(sig_d)

            if ipr_acc and ipr_acc != "-":
                ipr_ids[pid].add(ipr_acc)
            if ipr_d and ipr_d != "-":
                ipr_desc[pid].add(ipr_d)

            if go and go != "-":
                for term in go.split("|"):
                    term = term.strip()
                    if term:
                        go_terms[pid].add(term)

    return pfam_ids, ipr_ids, go_terms, sig_desc, ipr_desc

def join_set(s):
    return "-" if not s else ";".join(sorted(s))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--proteins_fa", required=True)
    ap.add_argument("--interpro_tsv", required=True)
    ap.add_argument("--swissprot_tsv", required=True)
    ap.add_argument("--uniref90_tsv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tx2gene = read_gene_transcript_map(args.gtf)
    prot_ids = read_fasta_ids(args.proteins_fa)

    swiss  = read_diamond_besthit(args.swissprot_tsv)
    uniref = read_diamond_besthit(args.uniref90_tsv)

    pfam_ids, ipr_ids, go_terms, sig_desc, ipr_desc = read_interproscan_tsv(args.interpro_tsv)

    header = [
        "gene_id",
        "transcript_id",
        "swissprot_sseqid","swissprot_evalue","swissprot_bitscore","swissprot_title",
        "uniref90_sseqid","uniref90_evalue","uniref90_bitscore","uniref90_title",
        "pfam_ids","interpro_ids","GO_terms",
        "signature_descriptions","interpro_descriptions"
    ]

    with open(args.out, "w") as out:
        out.write("\t".join(header) + "\n")
        for tid in prot_ids:
            gid = tx2gene.get(tid, "-")

            sp = swiss.get(tid, ("-","-","-","-"))
            ur = uniref.get(tid, ("-","-","-","-"))

            row = [
                gid, tid,
                sp[0], sp[1], sp[2], sp[3],
                ur[0], ur[1], ur[2], ur[3],
                join_set(pfam_ids.get(tid, set())),
                join_set(ipr_ids.get(tid, set())),
                join_set(go_terms.get(tid, set())),
                join_set(sig_desc.get(tid, set())),
                join_set(ipr_desc.get(tid, set())),
            ]
            out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    main()
