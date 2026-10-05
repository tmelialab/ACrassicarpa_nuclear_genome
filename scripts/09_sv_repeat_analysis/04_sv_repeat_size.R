library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
del <- read.delim(args[1]); ins <- read.delim(args[2])
outdir <- "output"; dir.create(outdir, showWarnings=FALSE)

del$SVTYPE <- "DEL"; ins$SVTYPE <- "INS"
sv <- rbind(del, ins)

sv$SIZE_CLASS <- cut(sv$SVLEN, c(50,100,500,1000,5000,10000,50000,Inf), right=FALSE,
 labels=c("50–99 bp","100–499 bp","500–999 bp","1–4.9 kb","5–9.9 kb","10–49.9 kb","≥50 kb"))

keep <- c("Unclassified","LTR_Gypsy","LTR_Copia","DNA_transposon","LINE","Simple_repeat","Non_repeat")
sv$CLASS <- ifelse(sv$DOMINANT_CLASS %in% keep, sv$DOMINANT_CLASS, "Other")

tab <- as.data.frame(table(sv$SVTYPE, sv$SIZE_CLASS, sv$CLASS))
names(tab) <- c("SVTYPE","SIZE_CLASS","CLASS","COUNT")
tab$PERCENT <- ave(tab$COUNT, tab$SVTYPE, tab$SIZE_CLASS, FUN=function(x) 100*x/sum(x))

write.table(tab, file.path(outdir,"sv_repeat_by_size.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cols <- c("Unclassified"="#999999","LTR_Gypsy"="#7B2CBF","LTR_Copia"="#E76F51",
          "DNA_transposon"="#E9C46A","LINE"="#2A9D8F","Simple_repeat"="#457B9D",
          "Non_repeat"="#D9D9D9","Other"="#8AB17D")

p <- ggplot(tab, aes(SIZE_CLASS, PERCENT, fill=CLASS)) +
  geom_col(width=0.8) +
  facet_wrap(~SVTYPE, ncol=1) +
  scale_fill_manual(values=cols) +
  labs(x="Structural variant size", y="Proportion of SVs (%)", fill="Dominant repeat class") +
  theme_classic(base_size=16) +
  theme(axis.text.x=element_text(angle=45,hjust=1), axis.title=element_text(face="bold"),
        legend.position="right", strip.background=element_blank(), strip.text=element_text(face="bold"))

ggsave(file.path(outdir,"sv_repeat_by_size.png"), p, width=8, height=6, dpi=600, bg="white")

# Median SV size for each repeat class
med <- aggregate(SVLEN ~ SVTYPE + CLASS, sv, function(x) c(N=length(x),Median=median(x),Mean=mean(x)))
med <- data.frame(SVTYPE=med$SVTYPE, CLASS=med$CLASS, N=med$SVLEN[,1],
                  MEDIAN_BP=round(med$SVLEN[,2]), MEAN_BP=round(med$SVLEN[,3]))
write.table(med, file.path(outdir,"sv_repeat_size_statistics.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

print(med)