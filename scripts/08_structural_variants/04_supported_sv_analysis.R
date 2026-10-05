# supported_sv_analysis.R

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
vcf <- args[1]; prefix <- args[2]

x <- readLines(vcf)
x <- x[!startsWith(x, "#")]
cat("VCF records:", length(x), "\n")

f <- strsplit(x, "\t", fixed = TRUE)
chrom <- vapply(f, `[`, "", 1)
pos <- as.numeric(vapply(f, `[`, "", 2))
info <- vapply(f, `[`, "", 8)

get_info <- function(x, key) {
  sapply(strsplit(x, ";", fixed = TRUE), function(z) {
    y <- z[startsWith(z, paste0(key, "="))]
    if (length(y)) sub(paste0("^", key, "="), "", y[1]) else NA
  })
}

sv <- data.frame(CHROM = chrom, POS = pos, SVTYPE = get_info(info, "SVTYPE"), SVLEN = abs(as.numeric(get_info(info, "SVLEN"))))
sv <- subset(sv, SVTYPE %in% c("INS", "DEL") & !is.na(SVLEN) & SVLEN >= 50)

sv$SIZE_CLASS <- cut(sv$SVLEN, breaks = c(50, 100, 500, 1000, 5000, 10000, 50000, Inf), right = FALSE, labels = c("50–99 bp", "100–499 bp", "500–999 bp", "1–4.9 kb", "5–9.9 kb", "10–49.9 kb", "≥50 kb"))

tab <- as.data.frame(table(sv$SIZE_CLASS, sv$SVTYPE))
names(tab) <- c("Size_class", "SV_type", "Count")
write.table(tab, paste0(prefix, "_size_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cols <- c("INS" = "#2E86AB", "DEL" = "#E76F51")

p1 <- ggplot(sv, aes(x = SVLEN, fill = SVTYPE)) +
  geom_histogram(bins = 80, alpha = 0.85) +
  scale_x_log10() +
  scale_fill_manual(values = cols) +
  facet_wrap(~SVTYPE, ncol = 1, scales = "free_y") +
  labs(x = "Structural variant size (bp, log scale)", y = "Number of variants") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", strip.background = element_rect(fill = "#F2F2F2", colour = NA), strip.text = element_text(face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(colour = "black"))

ggsave(paste0(prefix, "_size_distribution.png"), p1, width = 7, height = 6, dpi = 600, bg = "white")

p2 <- ggplot(tab, aes(x = Size_class, y = Count, fill = SV_type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = cols) +
  labs(x = "Structural variant size", y = "Number of variants", fill = "SV type") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top")
ggsave(paste0(prefix, "_size_classes.png"), p2, width = 8, height = 4, dpi = 600, bg = "white")

cat("\nINS/DEL variants:", nrow(sv), "\n")
print(tab)