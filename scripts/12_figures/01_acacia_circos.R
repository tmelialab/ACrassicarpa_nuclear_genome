rm(list = ls())

library(circlize)
library(Biostrings)
library(data.table)
library(IRanges)

# Usage:
# Rscript 01_acacia_circos.R \
#   <genome.fa> \
#   <braker.gtf> \
#   <repeatmasker.gff> \
#   <repeat_library.classified> \
#   <output_prefix>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    stop(
        paste(
            "Usage: Rscript 01_acacia_circos.R",
            "<genome.fa>",
            "<braker.gtf>",
            "<repeatmasker.gff>",
            "<repeat_library.classified>",
            "<output_prefix>"
        )
    )
}

genome_fasta   <- args[1]
gff_file       <- args[2]
repeat_file    <- args[3]
classified_file <- args[4]
output_prefix  <- args[5]

dir.create(
    dirname(output_prefix),
    recursive = TRUE,
    showWarnings = FALSE
)

top_n <- 25
window_size <- 1000000
top_n <- 25
window_size <- 1000000

# ============================================================
# Read genome and select largest contigs
# ============================================================

genome <- readDNAStringSet(genome_fasta)
names(genome) <- sub("\\s+.*$","",names(genome))

contigs <- data.frame(
  contig=names(genome),
  length=width(genome)
)

contigs <- contigs[order(contigs$length,decreasing=TRUE),]
top_contigs <- head(contigs,top_n)
selected_contigs <- top_contigs$contig

cat("Genome size:",sum(width(genome))/1e6,"Mb\n")
cat("Total contigs:",length(genome),"\n")
cat("Contigs plotted:",length(selected_contigs),"\n")
cat("Sequence represented:",sum(top_contigs$length)/1e6,"Mb\n")

# ============================================================
# Create genomic windows
# ============================================================

windows <- rbindlist(lapply(seq_len(nrow(top_contigs)),function(i){
  starts <- seq(1,top_contigs$length[i],by=window_size)
  data.table(
    chr=top_contigs$contig[i],
    start=starts,
    end=pmin(starts+window_size-1,top_contigs$length[i])
  )
}))

# ============================================================
# Read BRAKER annotation
# ============================================================

gff <- fread(
  gff_file,
  sep="\t",
  header=FALSE,
  quote="",
  fill=TRUE
)

gff <- gff[!grepl("^#",V1)]
gff <- gff[ncol(gff)>=9]

setnames(
  gff[1:9],
  c("chr","source","type","start","end","score","strand","phase","attributes")
)

gff <- gff[,1:9]
setnames(
  gff,
  c("chr","source","type","start","end","score","strand","phase","attributes")
)

gff[,start:=as.numeric(start)]
gff[,end:=as.numeric(end)]
gff <- gff[chr %in% selected_contigs]

# Prefer explicit gene features
genes <- gff[type=="gene",.(chr,start,end)]

# Fallback for GTF files without gene rows
if(nrow(genes)==0){
  tx <- gff[type %in% c("transcript","mRNA")]
  tx[,gene_id:=sub('.*gene_id[ ="]+([^";]+).*',"\\1",attributes)]
  
  if(nrow(tx)==0)
    stop("No gene or transcript features found in BRAKER annotation.")
  
  genes <- tx[,.(start=min(start),end=max(end)),by=.(chr,gene_id)]
  genes <- genes[,.(chr,start,end)]
}

genes[,midpoint:=(start+end)/2]

cat("Genes on plotted contigs:",nrow(genes),"\n")

# ============================================================
# Gene density
# ============================================================

windows[,gene_count:=0L]

for(i in seq_len(nrow(windows))){
  windows$gene_count[i] <- sum(
    genes$chr==windows$chr[i] &
      genes$midpoint>=windows$start[i] &
      genes$midpoint<=windows$end[i]
  )
}

# ============================================================
# Read RepeatMasker GFF
# ============================================================
rep <- read.delim(
  repeat_file,
  header=FALSE,
  sep="\t",
  comment.char="#",
  quote="",
  fill=TRUE,
  stringsAsFactors=FALSE
)

rep <- as.data.table(rep)

cat("Number of GFF columns:",ncol(rep),"\n")

if(ncol(rep)<9)
  stop("RepeatMasker GFF parsing failed: expected 9 columns.")

rep <- rep[,1:9]

setnames(
  rep,
  c(
    "chr","source","type","start","end",
    "score","strand","phase","attributes"
  )
)

rep[,start:=as.numeric(start)]
rep[,end:=as.numeric(end)]

rep <- rep[
  chr %in% selected_contigs &
    !is.na(start) &
    !is.na(end)
]

# Ensure start <= end
rep[,`:=`(
  start2=pmin(start,end),
  end2=pmax(start,end)
)]

rep[,`:=`(
  start=start2,
  end=end2
)]

rep[,c("start2","end2"):=NULL]

cat("Repeat records on plotted contigs:",nrow(rep),"\n")

head(rep)

# ============================================================
# MAP REPEATMODELER FAMILY IDs TO REPEAT CLASSES
# ============================================================



headers <- readLines(classified_file, warn=FALSE)
headers <- headers[grepl("^>",headers)]
headers <- sub("^>","",headers)

repeat_class <- data.table(
  family=sub("#.*$","",headers),
  class=sub("^[^#]*#","",headers)
)

# If a header has no # classification, mark it unknown
repeat_class[family==class,class:="Unknown"]

cat("Repeat families in classified library:",nrow(repeat_class),"\n")
print(head(repeat_class,10))

# Extract repeat-family ID from RepeatMasker GFF
rep[,family:=sub(
  '.*Motif:([^" ]+).*',
  '\\1',
  attributes
)]

# Assign class using consensi.fa.classified
rep[,repeat_class:=repeat_class$class[
  match(family,repeat_class$family)
]]

rep[is.na(repeat_class),repeat_class:="Unclassified"]

cat("\nRepeat classes found in GFF:\n")
print(sort(table(rep$repeat_class),decreasing=TRUE)[1:20])

# ============================================================
# GYPSY AND COPIA
# ============================================================

gypsy <- rep[
  grepl("Gypsy",repeat_class,ignore.case=TRUE)
]

copia <- rep[
  grepl("Copia",repeat_class,ignore.case=TRUE)
]

cat("\nGypsy records:",nrow(gypsy),"\n")
cat("Copia records:",nrow(copia),"\n")

# ============================================================
# Function to calculate nonredundant repeat coverage
# ============================================================

calculate_coverage <- function(repeat_table){
  
  result <- numeric(nrow(windows))
  
  for(i in seq_len(nrow(windows))){
    
    tmp <- repeat_table[
      chr==windows$chr[i] &
        start<=windows$end[i] &
        end>=windows$start[i]
    ]
    
    if(nrow(tmp)>0){
      
      clipped_start <- pmax(tmp$start,windows$start[i])
      clipped_end <- pmin(tmp$end,windows$end[i])
      
      rr <- reduce(
        IRanges(
          start=clipped_start,
          end=clipped_end
        )
      )
      
      result[i] <- 100*sum(width(rr))/
        (windows$end[i]-windows$start[i]+1)
    }
  }
  
  result
}

# ============================================================
# Repeat coverage tracks
# ============================================================

windows[,repeat_percent:=calculate_coverage(rep)]
windows[,gypsy_percent:=calculate_coverage(gypsy)]
windows[,copia_percent:=calculate_coverage(copia)]

# ============================================================
# GC content
# ============================================================

windows[,GC:=NA_real_]

for(i in seq_len(nrow(windows))){
  
  s <- subseq(
    genome[[windows$chr[i]]],
    start=windows$start[i],
    end=windows$end[i]
  )
  
  f <- alphabetFrequency(s,baseOnly=TRUE)
  denom <- sum(f[c("A","C","G","T")])
  
  if(denom>0)
    windows$GC[i] <- 100*(f["G"]+f["C"])/denom
}

cat("GC range:",range(windows$GC,na.rm=TRUE),"\n")
cat("Repeat coverage range:",range(windows$repeat_percent,na.rm=TRUE),"\n")
cat("Gypsy coverage range:",range(windows$gypsy_percent,na.rm=TRUE),"\n")
cat("Copia coverage range:",range(windows$copia_percent,na.rm=TRUE),"\n")

# ============================================================
# Save window statistics
# ============================================================

fwrite(
  windows,
  paste0(output_prefix,"_windows.tsv"),
  sep="\t"
)

# ============================================================
# Plot function
# ============================================================

draw_circos <- function(){
  
  circos.clear()
  
  circos.par(
    start.degree=90,
    gap.degree=4,
    track.margin=c(0.008,0.008),
    cell.padding=c(0,0,0,0),
    points.overflow.warning=FALSE
  )
  
 
  circos.initialize(
    factors=top_contigs$contig,
    xlim=cbind(
      rep(0,nrow(top_contigs)),
      top_contigs$length
    )
  )
  
  # ==========================================================
  # Track 1: contigs
  # ==========================================================
  
  circos.trackPlotRegion(
    ylim=c(0,1),
    track.height=0.085,
    bg.col="#E6E6E6",
    bg.border="white",
    panel.fun=function(x,y){
      
      sector <- CELL_META$sector.index
      len <- CELL_META$xlim[2]
      
      circos.axis(
        h="top",
        major.at=c(5,len),
        labels=c("0",sprintf("%.1f",len/1e6)),
        labels.cex=1,
        labels.niceFacing=TRUE,
        minor.ticks=0
      )
      
      # ctg_1_1 -> ctg1
      # ctg_48_1 -> ctg48
      label <- sub("^ctg_([0-9]+)_.*$", "ctg\\1", sector)
      
      circos.text(
        CELL_META$xcenter,
        0.5,
        label,
        facing="bending.inside",
        niceFacing=TRUE,
        cex=1,
        font=2
      )
    }
  )
  
  # ==========================================================
  # Track 2: gene density
  # ==========================================================
  gene_max <- max(windows$gene_count, na.rm=TRUE)
  if(gene_max == 0) gene_max <- 1
  
  circos.trackPlotRegion(
    ylim=c(0, gene_max),
    track.height=0.10,
    bg.col="#F7F7F7",
    bg.border="white",
    panel.fun=function(x, y){
      
      sector <- CELL_META$sector.index
      tmp <- windows[chr == sector]
      
      for(i in seq_len(nrow(tmp))){
        circos.rect(
          tmp$start[i], 0,
          tmp$end[i], tmp$gene_count[i],
          col="#3A923A",
          border=NA
        )
      }
      
      # Show vertical scale only on the first contig
      if(sector == top_contigs$contig[1]){
        circos.yaxis(
          side="left",
          at=c(0, gene_max),
          labels=c(0, gene_max),
          labels.cex=0.7,
          tick.length=0.02
        )
      }
    }
  )

  
  # ==========================================================
  # Track 3: total repeat coverage
  # ==========================================================
  
  circos.trackPlotRegion(
    ylim=c(0,100),
    track.height=0.10,
    bg.col="#F7F7F7",
    bg.border="white",
    panel.fun=function(x,y){
      
      sector <- CELL_META$sector.index
      tmp <- windows[chr==sector]
      
      for(i in seq_len(nrow(tmp))){
        circos.rect(
          tmp$start[i],0,
          tmp$end[i],tmp$repeat_percent[i],
          col="#C44E52",
          border=NA
        )
      }
      # Show vertical scale only on the first contig
      if(sector == top_contigs$contig[1]){
        circos.yaxis(
          side="left",
          at=c(5, 95),
          labels=c(0, 100),
          labels.cex=0.7,
          tick.length=0.02
        )
      }
    }
  )
  
  # ==========================================================
  # Track 4: GC content
  # ==========================================================
  
  gc_min <- floor(min(windows$GC,na.rm=TRUE))
  gc_max <- ceiling(max(windows$GC,na.rm=TRUE))
  
  if(gc_min==gc_max){
    gc_min <- gc_min-1
    gc_max <- gc_max+1
  }
  
  circos.trackPlotRegion(
    ylim=c(gc_min,gc_max),
    track.height=0.10,
    bg.col="#F7F7F7",
    bg.border="white",
    panel.fun=function(x,y){
      
      sector <- CELL_META$sector.index
      tmp <- windows[chr==sector]
      
      circos.lines(
        (tmp$start+tmp$end)/2,
        tmp$GC,
        col="#4C78A8",
        lwd=1.4
      )
      if(sector == top_contigs$contig[1]){
        circos.yaxis(
          side="left",
          at=c(gc_min+2, gc_max),
          labels=c(gc_min, gc_max),
          labels.cex=0.7,
          tick.length=0.02
        )
      }
    }
  )
  
  # ==========================================================
  # Track 5: Gypsy coverage
  # ==========================================================
  
  gypsy_max <- max(windows$gypsy_percent,na.rm=TRUE)
  if(!is.finite(gypsy_max) || gypsy_max==0) gypsy_max <- 1
  
  circos.trackPlotRegion(
    ylim=c(0,gypsy_max),
    track.height=0.085,
    bg.col="#F7F7F7",
    bg.border="white",
    panel.fun=function(x,y){
      
      sector <- CELL_META$sector.index
      tmp <- windows[chr==sector]
      
      for(i in seq_len(nrow(tmp))){
        circos.rect(
          tmp$start[i],0,
          tmp$end[i],tmp$gypsy_percent[i],
          col="#E89C5B",
          border=NA
        )
      }
      if(sector == top_contigs$contig[1]){
        circos.yaxis(
          side="left",
          at=c(5, round(gypsy_max)),
          labels=c(0, round(gypsy_max)),
          labels.cex=0.7,
          tick.length=0.02
        )
      }
    }
  )
  
  # ==========================================================
  # Track 6: Copia coverage
  # ==========================================================
  
  copia_max <- max(windows$copia_percent,na.rm=TRUE)
  if(!is.finite(copia_max) || copia_max==0) copia_max <- 1
  
  circos.trackPlotRegion(
    ylim=c(0,copia_max),
    track.height=0.085,
    bg.col="#F7F7F7",
    bg.border="white",
    panel.fun=function(x,y){
      
      sector <- CELL_META$sector.index
      tmp <- windows[chr==sector]
      
      for(i in seq_len(nrow(tmp))){
        circos.rect(
          tmp$start[i],0,
          tmp$end[i],tmp$copia_percent[i],
          col="#8E4585",
          border=NA
        )
      }
      if(sector == top_contigs$contig[1]){
        circos.yaxis(
          side="left",
          at=c(0, round(copia_max)),
          labels=c(0, round(copia_max)),
          labels.cex=0.7,
          tick.length=0.02
        )
      }
    }
  )

  # ==========================================================
  # Center
  # ==========================================================
  
  text(
    0,
    0.1,
    expression(italic("Acacia crassicarpa")),
    cex=2,
    font=2
  )
  
  # text(
  #   0,
  #   -0.02,
  #   "v4 genome",
  #   cex=0.9
  # )

  legend(
    -0.2,-0.03,
    inset=c(0,0),
    legend=c(
      "Gene density",
      "Repeat coverage",
      "GC content",
      "Gypsy",
      "Copia"
    ),
    fill=c(
      "#3A923A",
      "#C44E52",
      "#4C78A8",
      "#E89C5B",
      "#8E4585"
      
    ),
    border=NA,
    bty="n",
    cex=1,
    y.intersp=0.85
  )
  
  circos.clear()
}

# ============================================================
# PNG
# ============================================================

png(
  paste0(output_prefix,".png"),
  width=3000,
  height=3000,
  res=300
)

draw_circos()
dev.off()

# ============================================================
# PDF
# ============================================================

pdf(
  paste0(output_prefix,".pdf"),
  width=10,
  height=10,
  useDingbats=FALSE
)

draw_circos()
dev.off()

cat("\nSaved:\n")
cat(paste0(output_prefix,".png"),"\n")
cat(paste0(output_prefix,".pdf"),"\n")
cat(paste0(output_prefix,"_windows.tsv"),"\n")

