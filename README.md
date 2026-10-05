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

## Analysis workflows

The `scripts/` directory contains the computational workflows used for
assembly, annotation, comparative genomics, structural-variant analysis,
and gene-content analysis of the AcraUNRI *Acacia crassicarpa* genome.

Most compute-intensive workflows are provided as SLURM submission scripts.
Cluster-specific filesystem paths have been removed from the public versions.
Software modules or Conda environments may need to be adjusted for the
target HPC system.

### 01_raw_read_qc
Quality control of Oxford Nanopore long reads using NanoPlot.

### 02_genome_assembly
Long-read genome assembly using NextDenovo and iterative polishing using
NextPolish.

### 03_assembly_qc_purge
Assembly evaluation and redundancy reduction using BUSCO, Jellyfish,
Meryl, GenomeScope, Merqury, and purge_dups.

### 04_repeat_annotation
Species-specific repeat-library construction using RepeatModeler and repeat
annotation using RepeatMasker. The same species-specific repeat library was
applied to the Acra3RX reference genome for comparative repeat/SV analyses.

### 05_gene_annotation
RNA-seq-assisted BRAKER3 gene annotation and BUSCO assessment of the
predicted protein set.

### 06_functional_annotation
Functional annotation using InterProScan, Pfam, Swiss-Prot, and UniRef90.

### 07_genome_comparison
Whole-genome alignment of AcraUNRI against Acra3RX using minimap2.

The principal alignment settings were:

    minimap2 -x asm5 -c --eqx --secondary=no

Alignments of at least 50 kb were retained for the reported whole-genome
comparison.

### 08_structural_variants
Structural-variant discovery and read-level support analysis, including:

- SVIM-asm discovery
- ONT alignment against Acra3RX
- Sniffles2 force genotyping
- supported-SV extraction
- SV size and genotype analysis

### 09_sv_repeat_analysis
Analysis of repeat association with structural variants and scaffold- and
length-matched permutation analysis of deletion repeat enrichment.

### 10_sv_gene_analysis
Analysis of structural-variant overlap with genes and coding sequences,
including permutation-based tests for gene and CDS depletion.

### 11_orthology_pav
Orthogroup analysis and reciprocal genome-level validation of candidate
gene-content differences using OrthoFinder and miniprot.

### 12_figures
Scripts used to generate genome-level figures, including the AcraUNRI
circos representation and AcraUNRI–Acra3RX synteny visualization.

### 99_utilities
Utility scripts used during data preparation.

### Reproducibility

The scripts preserve the analysis parameters used in the manuscript.
Input and output paths have been generalized so that the workflows can be
run outside the original computing environment.
