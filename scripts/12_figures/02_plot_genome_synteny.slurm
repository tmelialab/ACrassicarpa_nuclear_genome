#!/bin/bash
#SBATCH --job-name=acra_asynt
#SBATCH --output=acra_asynt.%j.out
#SBATCH --error=acra_asynt.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00


set -euo pipefail

############################################
# USER SETTINGS
############################################

# Usage:
# sbatch 02_plot_genome_synteny.slurm \
#   <alignment.paf> \
#   <Acra3RX.fa> \
#   <AcraUNRI.fa> \
#   <output_dir> \
#   <asynt_directory>

PAF="${1:?Missing minimap2 PAF file}"
REF="${2:?Missing Acra3RX reference genome}"
QUERY="${3:?Missing AcraUNRI genome}"
OUTDIR="${4:?Missing output directory}"
ASYNT_DIR="${5:?Missing asynt directory}"

############################################
# VISUALIZATION PARAMETERS
############################################

# Do not use tiny alignments in the overview
MIN_ALIGNMENT=50000

# Adjacent same-orientation blocks may be joined
MAX_GAP=250000

# Final displayed block must span at least this much
MIN_BLOCK=250000

# Minimum unique sub-block retained during merging
MIN_SUBBLOCK=50000

############################################
# PREP
############################################

mkdir -p "$OUTDIR"

echo "=========================================="
echo "Acra3RX vs AcraUNRI alignment visualization"
echo "Started: $(date)"
echo "=========================================="

echo "PAF   : $PAF"
echo "REF   : $REF"
echo "QUERY : $QUERY"
echo "OUT   : $OUTDIR"

command -v samtools >/dev/null 2>&1 || {
    echo "ERROR: samtools not found." >&2
    exit 1
}

command -v Rscript >/dev/null 2>&1 || {
    echo "ERROR: Rscript not found." >&2
    exit 1
}

############################################
# CHECK FILES
############################################

for f in "$PAF" "$REF" "$QUERY" "$ASYNT_DIR/asynt.R"
do
    if [ ! -f "$f" ]; then
        echo "ERROR: file not found:"
        echo "$f"
        exit 1
    fi
done

############################################
# LOAD SOFTWARE
############################################

module load bioinformatics/samtools 2>/dev/null || true

############################################
# CREATE FASTA INDEXES
############################################

if [ ! -f "${REF}.fai" ]; then
    echo "Creating reference FASTA index"
    samtools faidx "$REF"
fi

if [ ! -f "${QUERY}.fai" ]; then
    echo "Creating query FASTA index"
    samtools faidx "$QUERY"
fi

############################################
# EXPORT SETTINGS TO R
############################################

export PAF REF QUERY ASYNT_DIR OUTDIR
export MIN_ALIGNMENT MAX_GAP MIN_BLOCK MIN_SUBBLOCK

############################################
# RUN R / ASYNT
############################################

Rscript - <<'RSCRIPT'

paf_file <- Sys.getenv("PAF")
ref_fai  <- paste0(Sys.getenv("REF"), ".fai")
qry_fai  <- paste0(Sys.getenv("QUERY"), ".fai")
asynt_r  <- file.path(Sys.getenv("ASYNT_DIR"), "asynt.R")
outdir   <- Sys.getenv("OUTDIR")

min_alignment <- as.numeric(Sys.getenv("MIN_ALIGNMENT"))
max_gap       <- as.numeric(Sys.getenv("MAX_GAP"))
min_block     <- as.numeric(Sys.getenv("MIN_BLOCK"))
min_subblock  <- as.numeric(Sys.getenv("MIN_SUBBLOCK"))

############################################
# LOAD ASYNT
############################################

source(asynt_r)

############################################
# IMPORT MINIMAP2 PAF
############################################

cat("\nImporting minimap2 PAF...\n")

alignments <- import.paf(paf_file)

cat("Raw alignment records:",
    nrow(alignments), "\n")

############################################
# IMPORT GENOME LENGTHS
############################################

ref_data <- import.genome(
    fai_file = ref_fai
)

query_data <- import.genome(
    fai_file = qry_fai
)

############################################
# FILTER SHORT ALIGNMENTS
############################################

alignments <- subset(
    alignments,
    Rlen >= min_alignment &
    Qlen >= min_alignment
)

cat(
    "Alignments >=",
    min_alignment,
    "bp:",
    nrow(alignments),
    "\n"
)

############################################
# MERGE NEARBY SAME-ORIENTATION ALIGNMENTS
#
# This is ONLY for visualization.
# We are NOT classifying rearrangements.
############################################

cat("\nGenerating merged alignment blocks...\n")

blocks <- get.synteny.blocks.multi(
    alignments,
    max_gap = max_gap,
    min_block_size = min_block,
    min_subblock_size = min_subblock
)

cat(
    "Merged alignment blocks:",
    nrow(blocks),
    "\n"
)

############################################
# CALCULATE BLOCK SIZE + ORIENTATION
############################################

blocks$ref_span <-
    abs(blocks$Rend - blocks$Rstart) + 1

blocks$query_span <-
    abs(blocks$Qend - blocks$Qstart) + 1

blocks$orientation <-
    ifelse(
        sign(blocks$Rend - blocks$Rstart) ==
        sign(blocks$Qend - blocks$Qstart),
        "+",
        "-"
    )

############################################
# SAVE EXACT BLOCKS USED FOR FIGURE
############################################

write.table(
    blocks,
    file = file.path(
        outdir,
        "merged_alignment_blocks.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

############################################
# CLEANER LABELS
############################################

# Change displayed AcraUNRI names:
# ctg_1_1 -> ctg1
clean_ctg <- function(x) {
    sub(
        "^ctg_([0-9]+)_.*$",
        "ctg\\1",
        x
    )
}

############################################
# PLOT FUNCTION
############################################

make_plot <- function(device_file, type="pdf") {

    pdf(
            device_file,
            width = 8,
            height = 3,
            pointsize = 11,
            useDingbats = FALSE
    )

    par(
        mar = c(2, 3, 2, 3),
        xpd = NA
    )

    ########################################
    # ASYNT PARALLEL ALIGNMENT PLOT
    ########################################

    plot.alignments.multi(
        blocks,

        reference_lens = ref_data$seq_len,
        query_lens     = query_data$seq_len,

        # automatically order/orient query contigs
        # to maximise correspondence
        no_reverse = FALSE,
        no_reorder = FALSE,

        # only sequences represented by blocks
        only_show_aligned_seqs = TRUE,

        # Acra3RX on top
        reference_above = TRUE,

        # cleaner broad ribbons
        sigmoid = TRUE,
        semitransparent = TRUE,

        show_outline = FALSE,

        # forward / reverse orientation
        cols = c(
            "#4C78A8",
            "#E45756"
        ),

        colour_by = "orientation",

        edge_width = 0.20,
        chrom_width = 0.07,

        # scaffold/contig labels
        show_labels = FALSE,
        angle_labels = TRUE,
        labels_cex = 0.55,
        labels_angle = 45,

        # show genome bars
        show_contigs = TRUE,

        centre = TRUE
    )

    mtext(
        "Acra3RX",
        side = 3,
        line = 0.2,
        adj = 0,
        font = 2
    )

    mtext(
        "AcraUNRI",
        side = 1,
        line = 0.2,
        adj = 0,
        font = 2
    )

#    legend(
#        "topright",
#        legend = c(
#            "Same orientation",
#            "Reverse orientation"
#        ),
#        fill = c(
#            "#4C78A880",
#            "#E4575680"
#        ),
#        border = NA,
#        bty = "n",
#        cex = 0.5
#    )

    dev.off()
}

############################################
# CREATE VECTOR FIGURES
############################################

pdf_file <- file.path(
    outdir,
    "Figure2_Acra3RX_vs_AcraUNRI_alignment_blocks_small.pdf"
)

make_plot(
    pdf_file,
    "pdf"
)

############################################
# SUMMARY
############################################

cat("\nReference sequences plotted:",
    length(unique(blocks$reference)),
    "\n")

cat("AcraUNRI contigs plotted:",
    length(unique(blocks$query)),
    "\n")

cat("Forward-orientation blocks:",
    sum(blocks$orientation == "+"),
    "\n")

cat("Reverse-orientation blocks:",
    sum(blocks$orientation == "-"),
    "\n")

cat("\nCreated:\n")
cat(pdf_file, "\n")

RSCRIPT

echo
echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="

ls -lh "$OUTDIR"