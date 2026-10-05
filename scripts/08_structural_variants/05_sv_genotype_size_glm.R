# sv_genotype_size.R
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1) stop("Usage: Rscript sv_genotype_size.R <VCF>")
vcf <- args[1]
outdir <- "output"
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)

# Read VCF
x <- readLines(vcf)
x <- x[!startsWith(x,"#")]
f <- strsplit(x,"\t",fixed=TRUE)
info <- vapply(f,`[`,"",8)
format <- vapply(f,`[`,"",9)
sample <- vapply(f,`[`,"",10)

get_info <- function(x,key) sapply(strsplit(x,";",fixed=TRUE),function(z){
  y <- z[startsWith(z,paste0(key,"="))]
  if(length(y)) sub(paste0("^",key,"="),"",y[1]) else NA
})

get_gt <- function(fmt,smp) mapply(function(a,b){
  a <- strsplit(a,":",fixed=TRUE)[[1]]
  b <- strsplit(b,":",fixed=TRUE)[[1]]
  i <- match("GT",a)
  if(!is.na(i) && i<=length(b)) b[i] else NA
},fmt,smp)

# Extract and filter SVs
sv <- data.frame(SVTYPE=get_info(info,"SVTYPE"),SVLEN=abs(as.numeric(get_info(info,"SVLEN"))),
                 GT=get_gt(format,sample),stringsAsFactors=FALSE)
sv$GT <- gsub("\\|","/",sv$GT)
sv <- subset(sv,SVTYPE %in% c("INS","DEL") & GT %in% c("0/1","1/0","1/1") &
               !is.na(SVLEN) & SVLEN>=50)
sv$GENOTYPE <- ifelse(sv$GT=="1/1","Homozygous-alt","Heterozygous")
sv$SVTYPE <- factor(sv$SVTYPE,levels=c("DEL","INS"))

# Size classes
size_labels <- c("50–99 bp","100–499 bp","500–999 bp","1–4.9 kb","5–9.9 kb","10–49.9 kb","≥50 kb")
sv$SIZE_CLASS <- cut(sv$SVLEN,breaks=c(50,100,500,1000,5000,10000,50000,Inf),
                     right=FALSE,labels=size_labels)

# Observed genotype proportions
tab <- as.data.frame(table(sv$SVTYPE,sv$SIZE_CLASS,sv$GENOTYPE))
names(tab) <- c("SV_type","Size_class","Genotype","Count")

tab$Percent <- ave(tab$Count,tab$SV_type,tab$Size_class,FUN=function(z){
  if(sum(z)==0) rep(NA_real_,length(z)) else 100*z/sum(z)
})

tab$SV_type <- factor(tab$SV_type,levels=c("DEL","INS"))
tab$Size_class <- factor(tab$Size_class,levels=size_labels)

write.table(tab,file.path(outdir,"sv_genotype_by_size.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Logistic regression
sv$HET <- as.integer(sv$GENOTYPE=="Heterozygous")
m <- glm(HET ~ log10(SVLEN)*SVTYPE,data=sv,family=binomial)

# Prediction at actual median SV length within each type × size class
pred <- aggregate(SVLEN ~ SVTYPE + SIZE_CLASS,data=sv,FUN=median)
pred <- pred[!is.na(pred$SIZE_CLASS),]
pred$SVTYPE <- factor(pred$SVTYPE,levels=c("DEL","INS"))
pred$SIZE_CLASS <- factor(pred$SIZE_CLASS,levels=size_labels)

pr <- predict(m,newdata=pred,type="link",se.fit=TRUE)
pred$Probability <- plogis(pr$fit)*100
pred$Lower <- plogis(pr$fit-1.96*pr$se.fit)*100
pred$Upper <- plogis(pr$fit+1.96*pr$se.fit)*100

names(pred)[names(pred)=="SVTYPE"] <- "SV_type"
pred$SV_type <- factor(pred$SV_type,levels=c("DEL","INS"))
write.table(pred,file.path(outdir,"sv_genotype_glm_predictions.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

# Save model output
sink(file.path(outdir,"sv_genotype_size_model.txt"))
cat("Total INS/DEL:",nrow(sv),"\n\n")
cat("Genotype counts:\n"); print(with(sv,table(SVTYPE,GENOTYPE)))
cat("\nLogistic regression:\n"); print(summary(m))
cat("\nOdds ratios:\n"); print(exp(coef(m)))
cat("\nPredicted heterozygosity by size class:\n")
print(pred[,c("SV_type","SIZE_CLASS","SVLEN","Probability","Lower","Upper")])
sink()

# Figure 2D
fill_cols <- c(
  "Heterozygous"="#4DAF4A",
  "Homozygous-alt"="#984EA3"
)

p <- ggplot(tab,aes(x=Size_class,y=Percent,fill=Genotype))+
  geom_col(position=position_dodge(width=0.80),width=0.70,na.rm=TRUE) +
  geom_errorbar(data=pred,aes(x=SIZE_CLASS,ymin=Lower,ymax=Upper),
                inherit.aes=FALSE,width=0.10,linewidth=0.55,colour="black") +
  geom_line(data=pred,aes(x=SIZE_CLASS,y=Probability,group=SV_type,colour="GLM-pred. P(HET)"),
            inherit.aes=FALSE,linewidth=1.1) +
  geom_point(data=pred,aes(x=SIZE_CLASS,y=Probability,colour="GLM-pred. P(HET)"),
             inherit.aes=FALSE,size=2.7) +
  facet_wrap(~SV_type,ncol=1) +
  scale_fill_manual(values=fill_cols) +
  scale_colour_manual(values=c("GLM-pred. P(HET)"="black")) +
  scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),
                     labels=function(x) paste0(x,"%")) +
  labs(x="Structural variant size",
       y="Genotype proportion / predicted heterozygosity (%)",
       fill="Obs. genotype",colour=NULL) +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        #axis.title=element_text(face="bold"),
        legend.position="top",
        legend.box="horizontal",
        strip.text=element_text(face="bold"),
        strip.background=element_rect(fill="white",colour="black"))

ggsave(file.path(outdir,"Figure2D_genotype_size_GLM.png"),
       p,width=10,height=8,dpi=600,bg="white")

cat("\nModel coefficients:\n"); print(summary(m)$coefficients)
cat("\nOdds ratios:\n"); print(exp(coef(m)))
cat("\nPredicted probabilities:\n")
print(pred[,c("SV_type","SIZE_CLASS","SVLEN","Probability","Lower","Upper")])
cat("\nCreated:",file.path(outdir,"Figure2D_genotype_size_GLM.png"),"\n")