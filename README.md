# *Acacia crassicarpa* Nuclear Genome

Assembly, polishing, QC, repeat annotation, structural + functional annotation pipeline
and wrap-up results for the *Acacia crassicarpa* nuclear genome.

> Dipindahkan dari `tmelialab/HPC` (Sep 2026): repo HPC kini hanya berisi panduan
> generik Mahameru; semua yang berbau *Acacia crassicarpa* tinggal di sini.

## Final stats (round-5 polish, Aug 2026)

| Item | Value |
|---|---|
| Assembly | 860,812,502 bp, 397 contigs (`ctg_1..ctg_397`, longest-first), N50 7 Mb |
| BUSCO embryophyta_odb10 (n=1614, miniprot) | C 99.8% [S 89.7% D 10.1%], F 1, M 3, E 23.6% |
| RepeatMasker | 65.49% masked (563.7 Mb) |
| BRAKER3 3.0.6 | 33,942 proteins |
| InterProScan 5.75 | 86% hit, 1,251 unique GO |
| tRNAscan-SE 2.0.12 | 666 high-confidence (+175 pseudo) |
| DIAMOND | Swiss-Prot 81.2% / UniRef90 98.8% (e ≤ 1e-5) |

## Layout

```
script/    SLURM + Python pipeline 01–14 (+90 SyRI bonus), paths sanitized ($HOME, /mgpfs/scratch/$USER)
results/   wrap-up: 01_assembly/ 02_busco/ 03_repeat/ 04_annotation/ + md5sum.txt
```

## Reproduce

1. Sesuaikan variabel path di tiap skrip (`GENOME`, `RNA_DIR`, `IPR`, `DB`, `*.fofn`) dengan akun HPC masing-masing.
2. conda env: `training_qc`, `repeats`, `syri_env`; partisi BRIN `short` (≤24 jam, ≤64 core) / `medium-small` (≤3 hari), `--mem` ≤240G.
3. Urutan: 01 → 02 (atau loop `02_run_polish_rounds.sh N`) → 03 → 04 → 05 → 06 → 07/08 → 09 → 10 → 11 → 12 → 13 → 14.
4. Detail per langkah + gotcha (LD_LIBRARY_PATH NextPolish, trailing newline `.fofn`, flag STAR `--outSAMstrandField intronMotif`, strip `*` sebelum InterProScan): lihat `script/README.md`.

## Data availability

File >100 MB tidak disimpan di git (batas GitHub): assembly FASTA + masked FASTA (~820–840 MB),
`braker.aa`/`gff3`/`codingseq`, `braker.tsv` (1.5G), BAM/reads. Yang ada di sini: `braker.gtf` (52M),
`genome_annotation_final.tsv` (12M), ringkasan BUSCO/repeat/DIAMOND/tRNA/InterPro + checksum.
Data penuh tersedia via lab (HPC BRIN Mahameru) — hubungi maintainer untuk akses/salinan.

## Catatan versi

- Snapshot pipeline lama + hasil assembly Maret 2026 sebelumnya ada di `tmelialab/HPC` (`script/`, `output/`, pra-Sep 2026).
- Repo ini = **rerun final**: polish round-5 + re-anotasi penuh 19–27 Agu 2026 (lihat `results/`).
