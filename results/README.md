# Results — *Acacia crassicarpa* nuclear genome (round-5 polish, Aug 2026)

Wrap-up hasil final. Checksum: `md5sum.txt`.

## `01_assembly/`

- `rename_map.tsv` — old NextPolish header → `ctg_1..ctg_397` (sorted longest-first) + length + rank.
- Assembly FASTA **tidak** disimpan di git (>800 MB/file): `acacia_final.fasta` (860,812,502 bp, 397 ctg)
  dan `acacia_final.fasta.masked` — tersedia di HPC/scratch lab (lihat README root § Data availability).

## `02_busco/` (embryophyta_odb10, n=1614, miniprot, BUSCO 6.0.0)

- `C:99.8% [S:89.7% D:10.1%] F:0.1% M:0.2%, E:23.6%` — `short_summary...txt` / `.json`, `busco_figure.png`.
- Polish tidak menaikkan kelengkapan (draft sudah jenuh); D ~10% = haplotig/heterozigositas tertahan.

## `03_repeat/`

- `acacia_final.fasta.tbl` — RepeatMasker: **65.49% masked** (563.7 Mb), Unclassified 64.53% + Simple 0.95%, GC 35.62%.
  (Library = `consensi.fa` agregat round 1–5, belum RepeatClassifier → semua repeat "Unclassified".)

## `04_annotation/` (33,942 protein BRAKER3 3.0.6)

- `braker.gtf` (52M) — model gen final.
- `genome_annotation_final.tsv` (12M, 33,943 baris) — tabel master: koordinat GTF + Pfam/InterPro/GO + best-hit Swiss-Prot & UniRef90.
- `diamond_swissprot.tsv` (27,564 hit, 81.2%) / `diamond_uniref90.tsv` (33,528 hit, 98.8%), evalue ≤1e-5.
- `trna.out` / `trna.fa` / `trna.stats` — tRNAscan-SE 2.0.12: 841 kandidat, **666 high-confidence** (tanpa pseudo).
- `interproscan_summary.txt` / `pfam.txt` / `top20_pfam.txt` — 86% protein ada hit, 54% ber-GO, 1,251 unique GO.
- Protein FASTA penuh (`braker.aa` 15M) dan `braker.tsv` 1.5G **tidak** di git (besar) — di HPC lab.
