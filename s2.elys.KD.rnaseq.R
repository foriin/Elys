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

sal.ekd.d <- dir("/home/artem/IMG/Projects/Elys/S2_RNAseq/salmon_aln_pc/",
             full.names = T)
sal.ekd.f <- sal.ekd.f <- file.path(sal.ekd.d, "quant.genes.sf")

names(sal.ekd.f) <- paste0(rep(c("ZKD", "ElysKD"), each = 2), 2:3)

txi.s2.elys <- tximport(sal.ekd.f, "salmon", txOut = T, importer=read.delim,
                        countsFromAbundance = "scaledTPM")

# Make heatmaps of TPM correlations for each experiment
MainCorrelations(txi.s2.elys$counts %>% as.data.frame(),
                 corr.desc = "ELYS_reads", file.p = "plots")
MainCorrelations(txi.s2.lam$counts %>% as.data.frame(),
                 corr.desc = "LAM_reads", file.p = "plots")

# Merge TPM data with genes table



s2.elys.TPM <- merge(dm3.genes, txi.s2.elys$abundance %>% as.data.frame() %>%
                      mutate(id = rownames(txi.s2.elys$abundance)),
                     by = "id") %>% 
  mutate(ZTPMav = rowMeans(select(., matches("ZKD"))),
         ElysTPMav = rowMeans(select(., matches("ElysKD"))),
         chrom = ifelse(chr == "X", "X", "A")) %>% 
  select(-matches("KD"))

 
# DESEQ2

# Elys
samples <- data.frame(row.names = names(sal.ekd.f),
                      conditions = rep(c("wt", "elysKD"), each = 2))

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
minz <- min(s2.elys.TPM$ZTPMav[s2.elys.TPM$ZTPMav != 0])
res.df.elyskd <- as.data.frame(res.elyskd) %>% mutate(id = rownames(res.elyskd))
res.df.elyskd <- merge(s2.elys.TPM, res.df.elyskd, by = "id", all.x = T) %>% 
  select(1:10, 17) %>% mutate(log2FC = ifelse(ZTPMav != 0,
                                             log2(ElysTPMav/ZTPMav),
                                                  log2(ElysTPMav/minz)), 
                             diffex = ifelse(padj < 0.05, 1, 0))
# save(res.df.lamkd, res.df.elyskd, file = "RData/s2.lam.elys.deseq.RData")
write.xlsx(res.df.elyskd, "tables/s2.elysKD.difex.xlsx", overwrite = T)
save(s2.elys.TPM, res.df.elyskd, file = "RData/s2.elys.difex.tpm.RData")
save(txi.s2.elys, file = "RData/s2.elysKD.txi.RData")

rde <- res.df.elyskd %>% filter(padj < 0.05)
test <- makeGRangesFromDataFrame(rde, seqnames.field = "chr")
test$score <- rde$log2FoldChange
export.bedGraph(test, "bed/elyskd.foldchange.bedgraph")
