# Genome Assembly, Polishing, QC, Annotation & Comparative Pipeline — *Acacia crassicarpa*

Skrip bernomor sesuai urutan eksekusi. Semua path absolut asal sudah disanitasi:
`$HOME` = home HPC masing-masing, `/mgpfs/scratch/$USER` = scratch masing-masing.
Sesuaikan blok `USER SETTINGS` / variabel `GENOME`, `RNA_DIR`, `IPR`, `DB` sebelum submit.

## EN

1. **De novo assembly & pre-QC** — `01_assembly_nextdenovo.sh` (`01_nextdenovo.cfg`, `01_input.fofn`):
   NanoPlot QC on ONT reads, then `nextDenovo` → `hasil/03.ctg_graph/nd.asm.fasta`.
2. **Short-read polishing** — `02_polish_nextpolish.sh` (`02_polish.cfg`, `02_sgs.fofn`) single round;
   `02_run_polish_rounds.sh` loops N rounds (`task=best`, `rerun=5` = job retries, not polish iterations).
   Final assembly used here: round 5.
3. **QC** — `03_qc_busco.sh` (template) and `03_run_busco_polish.sbatch` (executed version):
   `busco -m genome -l embryophyta_odb10 --offline`, miniprot mode, + `busco --plot`.
4. **Repeat prep** — `04_repeatmodeler_builddb.sh`: `BuildDatabase -name acacia_db` (RepeatModeler 2.0.7, no `-engine` flag).
5. **Repeat discovery** — `05_run_repeatmodeler.sh`: `RepeatModeler -database acacia_db -threads 60 -LTRStruct` (`-pa` deprecated).
6. **Repeat masking** — `06_run_repeatmasker.sh`: soft-mask (`-xsmall`), auto-fallback to `RM_*/consensi.fa` if `acacia_db-families.fa` absent.
7. **RNA data prep** — `07_download_sra_reads.sh`: `fasterq-dump` list in `sra_dump.txt` (kept for reproducibility; rerun reused existing fastqs).
8. **RNA alignment** — `08_align_rnaseq_star.sh`: STAR index (`--genomeSAindexNbases 13`) + 4 samples → sorted BAM **with `--outSAMstrandField intronMotif`** (required for BRAKER3/StringTie; earlier run without it was redone).
9. **Gene prediction** — `09_run_braker3.sh`: BRAKER3 (`braker.pl` 3.0.6) with RNA BAMs + Viridiplantae proteins, `--softmasking`. Needs `GENEMARK_PATH`, `PROTHINT_PATH`, `~/.gm_key`.
10. **Domains & GO** — `10_run_interproscan.sh`: InterProScan 5.75 (`-appl Pfam,SMART,ProSiteProfiles,ProSitePatterns,SUPERFAMILY`), 64 CPU. Input must have `*` (stop) stripped → `braker.clean.aa`.
11. **tRNA** — `11_run_trnascan.sh`: `tRNAscan-SE -E` (Eukaryotic).
12. **Swiss-Prot homology** — `12_run_diamond_swissprot.sh`: `diamond blastp --sensitive -e 1e-5`, 12-col outfmt.
13. **UniRef90 homology** — `13_run_diamond_uniref90.sh`: same params, 86G DB (Query-indexed mode).
14. **Master table** — `14_make_genome_annotation_tsv.py`: merges GTF + InterPro + DIAMOND best hits → `genome_annotation_final.tsv`. Handles 12-col DIAMOND outfmt (evalue col 10, bitscore col 11).
90. **Comparative (bonus)** — `90_run_syri_sv_analysis.sh`: minimap2 `--eqx` → SyRI → plotsr.

## ID

(Versi Indonesia — ringkas; detail sama dengan EN di atas.)

1. **Perakitan** (`01_*`): NanoPlot lalu NextDenovo.
2. **Polish** (`02_*`): NextPolish 1 round atau loop N round via `02_run_polish_rounds.sh`.
3. **QC** (`03_*`): BUSCO lineage `embryophyta_odb10` offline + plot.
4–6. **Repeat** (`04–06`): BuildDatabase → RepeatModeler → RepeatMasker soft-mask.
7. **Unduh RNA** (`07`): `fasterq-dump` (arsip; rerun memakai fastq yang sudah ada).
8. **Align RNA** (`08`): STAR + flag wajib `--outSAMstrandField intronMotif`.
9. **BRAKER3** (`09`): RNA BAM + protein Viridiplantae; butuh `GENEMARK_PATH`, `PROTHINT_PATH`, `~/.gm_key`.
10. **InterProScan** (`10`): hapus `*` dari `braker.aa` dulu.
11. **tRNAscan** (`11`): mode Eukaryotic.
12–13. **DIAMOND** (`12–13`): Swiss-Prot & UniRef90, evalue 1e-5.
14. **Tabel master** (`14`): gabung GTF + InterPro + DIAMOND.
90. **Komparatif bonus** (`90`): minimap2 → SyRI → plotsr.

## Konvensi

- conda env: `training_qc` (assembly/polish/STAR/BRAKER/diamond), `repeats` (RepeatModeler/Masker), `syri_env` (SyRI).
- Partisi BRIN: `short` (≤24 jam, ≤64 core/job), `medium-small` (≤3 hari); `--mem` maks ~240G (256G ditolak scheduler).
- NextPolish butuh `LD_LIBRARY_PATH=$HOME/.conda/envs/training_qc/lib` **di dalam** job script.
- Setiap baris `.fofn` wajib diakhiri newline (bug `seq_split` memotong huruf terakhir bila tidak).
