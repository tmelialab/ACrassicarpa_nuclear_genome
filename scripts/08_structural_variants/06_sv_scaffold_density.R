# sv_scaffold_density.R
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]; fai <- args[2]; prefix <- args[3]

x <- readLines(vcf); x <- x[!startsWith(x, "#")]
f <- strsplit(x, "\t", fixed = TRUE)
chrom <- vapply(f, `[`, "", 1)
info <- vapply(f, `[`, "", 8)

get_info <- function(x, key) sapply(strsplit(x, ";", fixed = TRUE), function(z) {
  y <- z[startsWith(z, paste0(key, "="))]
  if(length(y)) sub(paste0("^", key, "="), "", y[1]) else NA
})

sv <- data.frame(CHROM = chrom, SVTYPE = get_info(info, "SVTYPE"))
sv <- subset(sv, SVTYPE %in% c("INS", "DEL"))

len <- read.table(fai, sep = "\t", header = FALSE)
len <- len[,1:2]; names(len) <- c("CHROM", "LENGTH")

count <- as.data.frame(table(sv$CHROM, sv$SVTYPE))
names(count) <- c("CHROM", "SVTYPE", "COUNT")

d <- merge(count, len, by = "CHROM")
d$MB <- d$LENGTH / 1e6
d$DENSITY <- d$COUNT / d$MB

# Focus figure on major scaffolds
major <- unique(d$CHROM[d$LENGTH >= 1e6])
d2 <- subset(d, CHROM %in% major)

# keep numerical scaffold order
d2$N <- as.numeric(sub("scaffold_", "", d2$CHROM))
d2 <- d2[order(d2$N),]
d2$CHROM <- factor(d2$CHROM, levels = unique(d2$CHROM))

cols <- c("INS" = "#2E86AB", "DEL" = "#E76F51")

p <- ggplot(d2, aes(CHROM, DENSITY, fill = SVTYPE)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols) +
  labs(x = "Acra3RX scaffold", y = "ONT-supported SVs/Mb", fill = "SV type") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top")

ggsave(paste0(prefix, "_sv_density.png"), p, width = 8, height = 4, dpi = 600, bg = "white")
write.table(d, paste0(prefix, "_sv_density.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

print(d2)