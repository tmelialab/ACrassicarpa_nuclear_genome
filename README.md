# AcraUNRI — Acacia crassicarpa nuclear genome

This repository contains the genome assembly, gene annotation, structural
variant dataset, and computational workflows associated with the AcraUNRI
*Acacia crassicarpa* genome study.

## Manuscript

**Long-Read Comparative Genomics Reveals Repeat-Associated Structural Variation and Conserved Gene Content in Acacia crassicarpa**

The study presents a long-read genome assembly of *Acacia crassicarpa* and
a comparative genomic analysis against the previously published Acra3RX
assembly.

## AcraUNRI genome assembly

| Metric | Value |
|---|---:|
| Assembly size | 831,473,243 bp |
| Contigs | 349 |
| Contig N50 | 7,801,588 bp |
| Longest contig | 23,730,391 bp |
| GC content | 35.64% |
| Complete BUSCO | 99.8% |
| Single-copy BUSCO | 91.7% |
| Duplicated BUSCO | 8.1% |
| Repeat content | 64.78% |
| Protein-coding genes | 29,789 |
| Predicted transcripts | 33,363 |
| Protein-set BUSCO | 98.5% |

## Repository contents

```text
.
├── genome/
│   ├── AcraUNRI_v1.softmasked.fa.gz
│   ├── AcraUNRI_v1.softmasked.fa.gz.fai
│   └── AcraUNRI_v1.softmasked.fa.gz.gzi
│
├── annotation/
│   ├── AcraUNRI_braker.gtf.gz
│   └── AcraUNRI_braker.aa.gz
│
├── supplementary/
│   ├── Supplementary_File_1_AcraUNRI_Acra3RX_ONT_supported_SVs.vcf.gz
│   └── Supplementary_File_1_AcraUNRI_Acra3RX_ONT_supported_SVs.vcf.gz.tbi
│
├── scripts/
│   ├── 01_raw_read_qc/
│   ├── 02_genome_assembly/
│   ├── 03_assembly_qc_purge/
│   ├── 04_repeat_annotation/
│   ├── 05_gene_annotation/
│   ├── 06_functional_annotation/
│   ├── 07_genome_comparison/
│   ├── 08_structural_variants/
│   ├── 09_sv_repeat_analysis/
│   ├── 10_sv_gene_analysis/
│   ├── 11_orthology_pav/
│   ├── 12_figures/
│   └── 99_utilities/
│
└── results/
