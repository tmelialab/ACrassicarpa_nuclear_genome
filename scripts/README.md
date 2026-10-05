# AcraUNRI analysis workflows

This directory contains the computational workflows used for assembly,
annotation, comparative genomics, structural-variant analysis, and
gene-content analysis of the AcraUNRI *Acacia crassicarpa* genome.

The scripts are organized according to the major analysis stages described
in the manuscript.

## Directory structure

### 01_raw_read_qc
Quality control of Oxford Nanopore long reads using NanoPlot.

### 02_genome_assembly
Long-read genome assembly using NextDenovo and iterative polishing using
NextPolish.

### 03_assembly_qc_purge
Assembly evaluation and redundancy reduction, including:

- BUSCO
- Jellyfish
- Meryl
- GenomeScope
- Merqury
- purge_dups

### 04_repeat_annotation
Species-specific repeat-library construction with RepeatModeler and repeat
annotation with RepeatMasker.

The same species-specific repeat library was also applied to the Acra3RX
reference genome for comparative repeat/SV analyses.

### 05_gene_annotation
RNA-seq-assisted BRAKER3 gene annotation and BUSCO evaluation of the
predicted protein set.

### 06_functional_annotation
Functional annotation using:

- InterProScan
- Pfam
- Swiss-Prot
- UniRef90

### 07_genome_comparison
Whole-genome alignment of AcraUNRI against Acra3RX using minimap2.

The principal alignment used:

    minimap2 -x asm5 -c --eqx --secondary=no

Alignments of at least 50 kb were retained for the reported whole-genome
comparison.

### 08_structural_variants
Structural-variant discovery and read-level support analysis.

The workflow includes:

- SVIM-asm discovery
- ONT alignment against Acra3RX
- Sniffles2 force genotyping
- supported-SV extraction
- SV size/genotype analysis

### 09_sv_repeat_analysis
Analysis of repeat association with structural variants and matched
permutation analysis of deletion repeat enrichment.

### 10_sv_gene_analysis
Analysis of structural-variant overlap with genes and coding sequences,
including matched permutation tests for gene/CDS depletion.

### 11_orthology_pav
Orthogroup analysis and reciprocal genome-level validation of candidate
gene-content differences using OrthoFinder and miniprot.

### 12_figures
Scripts used to generate genome-level figures, including the AcraUNRI
circos representation and AcraUNRI-Acra3RX synteny visualization.

### 99_utilities
Utility scripts used during data preparation.

## HPC environment

Most compute-intensive workflows are provided as SLURM submission scripts.
Cluster-specific filesystem paths have been removed from the public versions.

Software modules or Conda environments may need to be adjusted for the
target HPC system.

## Reproducibility note

The scripts preserve the analysis parameters used for the manuscript.
Input and output paths have been generalized so that the workflows can be
run outside the original computing environment.
