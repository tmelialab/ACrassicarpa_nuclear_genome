# sv_window_density.R
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]; fai <- args[2]
outdir <- "output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

win <- ifelse(length(args) >= 3, as.numeric(args[3]) * 1e6, 1e6)
min_scaffold <- 1e6

ref <- read.table(fai, sep="\t")[,1:2]; names(ref) <- c("CHROM","LENGTH")
ref <- ref[ref$LENGTH >= min_scaffold,]
ref$OFFSET <- c(0, head(cumsum(ref$LENGTH), -1))
ref$MIDGENOME <- (ref$OFFSET + ref$LENGTH/2) / 1e6

x <- readLines(vcf); x <- x[!startsWith(x, "#")]
f <- strsplit(x, "\t", fixed=TRUE)
sv <- data.frame(CHROM=vapply(f, `[`, "", 1), POS=as.numeric(vapply(f, `[`, "", 2)), INFO=vapply(f, `[`, "", 8))

get_info <- function(x, key) sapply(strsplit(x, ";", fixed=TRUE), function(z) { y <- z[startsWith(z, paste0(key,"="))]; if(length(y)) sub(paste0("^",key,"="), "", y[1]) else NA })

sv$SVTYPE <- get_info(sv$INFO, "SVTYPE")
sv <- subset(sv, SVTYPE %in% c("INS","DEL") & CHROM %in% ref$CHROM)
sv$WIN <- floor((sv$POS - 1) / win) + 1

w <- do.call(rbind, lapply(seq_len(nrow(ref)), function(i) {
  s <- seq(1, ref$LENGTH[i], by=win); e <- pmin(s+win-1, ref$LENGTH[i])
  data.frame(CHROM=ref$CHROM[i], WIN=seq_along(s), START=s, END=e, SIZE=e-s+1, OFFSET=ref$OFFSET[i])
}))

cnt <- aggregate(POS ~ CHROM + WIN + SVTYPE, sv, length); names(cnt)[4] <- "COUNT"
d <- rbind(transform(w, SVTYPE="DEL"), transform(w, SVTYPE="INS"))
d <- merge(d, cnt, by=c("CHROM","WIN","SVTYPE"), all.x=TRUE)
d$COUNT[is.na(d$COUNT)] <- 0
d$DENSITY <- d$COUNT / (d$SIZE/1e6)
d$GENOME_MB <- (d$OFFSET + (d$START+d$END)/2) / 1e6

write.table(d, file.path(outdir, "acra_supported_1Mb_windows.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cols <- c("DEL"="#E76F51", "INS"="#2E86AB")
p <- ggplot(d, aes(GENOME_MB, DENSITY, colour=SVTYPE, group=interaction(CHROM,SVTYPE))) +
  geom_line(linewidth=0.55) +
  facet_wrap(~SVTYPE, ncol=1, scales="free_y") +
  scale_colour_manual(values=cols) +
  scale_x_continuous(breaks=ref$MIDGENOME, labels=ref$CHROM, expand=c(0.005,0.005)) +
  labs(x="Acra3RX scaffold", y="ONT-supported SVs per Mb") +
  theme_classic(base_size=13) +
  theme(legend.position="none", strip.background=element_blank(), strip.text=element_text(face="bold"), axis.text.x=element_text(angle=60, hjust=1), axis.title=element_text(face="bold"))

ggsave(file.path(outdir, "acra_supported_1Mb_sv_density.png"), p, width=12, height=6, dpi=600, bg="white")

cat("Output folder:", normalizePath(outdir), "\n")
cat("Window size:", win/1e6, "Mb\n")
cat("Scaffolds:", nrow(ref), "\n")
cat("INS:", sum(sv$SVTYPE=="INS"), "\n")
cat("DEL:", sum(sv$SVTYPE=="DEL"), "\n")