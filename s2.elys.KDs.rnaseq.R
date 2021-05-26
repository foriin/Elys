library(openxlsx)
library(dm3)
library(GenomicRanges)
library(dplyr)
library(tximport)
library(ggplot2)
library(gplots)
library(DESeq2)

MakeCorMat <- function(data, cormethod, use.vals = "complete.obs"){
  cormat <- as.matrix(data) %>%
    cor(method = cormethod, use = use.vals) %>% 
    round(digits = 2)
}

MainCorrelations <- function(data, use.opt="everything", corr.desc, createPDF=T,
                             file.p, ...) {  
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


# Import output of Salmon
sal.lkd.d <- dir("/home/artem/IMG/Projects/LAM.SpG.SpC/s2_RNA-seq/LamKD/salmon_out/",
             full.names = T)
sal.lkd.f <- file.path(sal.lkd.d, "quant.genes.sf")

names(sal.lkd.f) <- c("ZKD1", "ZKD2", "LamKD1", "LamKD2")

txi.s2.lam <- tximport(sal.lkd.f, "salmon", txOut = T, importer=read.delim,
                       countsFromAbundance = "scaledTPM")

sal.ekd.d <- dir("/home/artem/IMG/Projects/LAM.SpG.SpC/s2_RNA-seq/ElysKD/salmon_out",
             full.names = T)
sal.ekd.f <- sal.ekd.f <- file.path(sal.ekd.d, "quant.genes.sf")

names(sal.ekd.f) <- paste0(rep(c("ZKD", "ElysKD"), each = 3), 1:3)

txi.s2.elys <- tximport(sal.ekd.f, "salmon", txOut = T, importer=read.delim,
                        countsFromAbundance = "scaledTPM")

# Make heatmaps of TPM correlations for each experiment
MainCorrelations(txi.s2.elys$counts %>% as.data.frame(),
                 corr.desc = "ELYS_reads", file.p = "plots")
MainCorrelations(txi.s2.lam$counts %>% as.data.frame(),
                 corr.desc = "LAM_reads", file.p = "plots")

# Merge TPM data with genes table

s2.lam.TPM <- merge(dm3.genes, txi.s2.lam$abundance %>% as.data.frame() %>%
  mutate(id = rownames(txi.s2.lam$abundance)), by = "id")%>% 
  mutate(WT_TPM = rowMeans(select(., matches("ZKD"))),
         LamKD_TPM = rowMeans(select(., matches("LamKD"))),
         chrom = ifelse(chr != "X", "A", "X")) %>% 
  select(-matches("KD[12]"))
s2.lam.TPM$chrom <- factor(s2.lam.TPM$chrom, levels = c("X", "A"))
# write.xlsx(s2.lam.TPM, file = "tables/s2.rnaseq.xlsx")


s2.elys.TPM <- merge(dm3.genes, txi.s2.elys$abundance %>% as.data.frame() %>%
                      mutate(id = rownames(txi.s2.elys$abundance)),
                     by = "id") %>% 
  mutate(ZTPMav = rowMeans(select(., matches("ZKD"))),
         ElysTPMav = rowMeans(select(., matches("ElysKD"))),
         chrom = ifelse(chr == "X", "X", "A")) %>% 
  select(-matches("KD"))

# save(s2.lam.TPM,s2.elys.TPM, file = "RData/s2.tpms.RData")


# DESEQ2
# Lam
samples <- data.frame(row.names = names(sal.lkd.f),
                      conditions = rep(c("wt", "lamKD"), each = 2))

ddsTxi.l <- DESeqDataSetFromTximport(txi.s2.lam, colData = samples,
                                   design =~ conditions)

keep <- rowSums(counts(ddsTxi.l)) >= 10

ddsTxi.l <- ddsTxi.l[keep,]
ddsTxi.l$conditions <- relevel(ddsTxi.l$conditions, ref = "wt")

dds.l <- DESeq(ddsTxi.l)
res.lamkd <- results(dds.l)
resultsNames(dds.l)

resLFC <- lfcShrink(dds.l, coef = "conditions_lamKD_vs_wt", type = "apeglm",
                    lfcThreshold = 1)
pdf("./plots/MA.plot.pdf")
DESeq2::plotMA(resLFC, alpha = 0.05, ylim=c(-6, 6))
dev.off()

summary(res.lamkd, alpha = 0.05)
res.df.lamkd <- as.data.frame(res.lamkd) %>% mutate(id = rownames(res.lamkd))
res.df.lamkd <- merge(dm3.genes, res.df.lamkd, by = "id", all.y = T)


# Elys
samples <- data.frame(row.names = names(sal.ekd.f),
                      conditions = rep(c("wt", "elysKD"), each = 3))

ddsTxi.e <- DESeqDataSetFromTximport(txi.s2.elys, colData = samples,
                                     design =~ conditions)

keep <- rowSums(counts(ddsTxi.e)) >= 10

ddsTxi.e <- ddsTxi.e[keep,]
ddsTxi.e$conditions <- relevel(ddsTxi.e$conditions, ref = "wt")

dds.e <- DESeq(ddsTxi.e)
res.elyskd <- results(dds.e)
resultsNames(dds.e)

resLFC <- lfcShrink(dds.e, coef = "conditions_elysKD_vs_wt", type = "apeglm",
                    lfcThreshold = 1)
pdf("./plots/MA.plot.pdf")
DESeq2::plotMA(resLFC, alpha = 0.05, ylim=c(-6, 6))
dev.off()

summary(res.elyskd, alpha = 0.05)
res.df.elyskd <- as.data.frame(res.elyskd) %>% mutate(id = rownames(res.elyskd))
res.df.elyskd <- merge(dm3.genes, res.df.elyskd, by = "id", all.y = T)
# save(res.df.lamkd, res.df.elyskd, file = "RData/s2.lam.elys.deseq.RData")

res.df.elyskd



