library(tximport)
library(dplyr)
library(GenomicRanges)
library(data.table)
library(gplots)

MakeCorMat <- function(data, cormethod, use.vals = "complete.obs"){
  cormat <- as.matrix(data) %>%
    cor(method = cormethod, use = use.vals) %>% 
    round(digits = 2)
}

MainCorrelations <- function(data, use.opt="everything", corr.desc, createPDF=T, file.p, ...) {  
  lapply(c("spearman", "pearson"), function(meth){
    cors <- MakeCorMat(data, meth)
    assign(paste0(corr.desc, ".", meth, ".cor"), cors)
    if (createPDF == T){
      options(warn=-1)
      pdf(file=file.path(file.p, paste0(corr.desc,  "_", meth, "_correlation_heatmap", ".pdf")), width=12, height=12)
      heatmap.2(cors, col=greenred(200), breaks=seq(from=-1, to=1, by=0.01),
                Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=cors,
                notecol="white", notecex = 1.5, margins=c(7,7),
                main=paste(meth, "'s correlation coefficients and hierarchical clustering for\n", "'",
                           sub("\\d+\\.(.*)", "\\1", corr.desc), "'", sep=""),
                cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2),
                lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="density")
      options(warn=0)
      dev.off()
    }  
  })
}


files <- paste0("/home/artem/IMG/Projects/Elys/RNA-seq_10.07.19/",
                "larva_testes_RNA-seq_Elys_Lam_KD/salmon_out/sample",
                1:9)
files <- file.path(files, "quant.genes.sf")
names(files) <- paste0(rep(c("WT", "ElysKD", "LamKD"), each = 3), 1:3)
txi.lartest <- tximport(files, "salmon", txOut = T, importer=read.delim)

cor(txi.lartest$abundance[,7], txi.lartest$abundance[,9])
MainCorrelations(txi.lartest$abundance %>% as.data.frame(), corr.desc = "TPMs", file.p = ".")

