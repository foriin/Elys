library(GenomicRanges)
library(DESeq2)
library(data.table)
library(dplyr)
library(stringr)
library(gplots)
library(tximport)
library(edgeR)
library(reshape2)
library(openxlsx)

genes <- fread("~/Database/dmel/Annotations/dmel-gene-r5.57.gff") %>% 
  dplyr::select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('^ID=(FBgn[0-9]+);.*', '\\1', attr),
         chr = paste0("chr", chr), gene_name = sub('.*Name=([^;]+);.*', '\\1', attr),
         tss = ifelse(strand == "+", start, end)) %>% 
  dplyr::select(-attr)

genes.gr <- GRanges(
  seqnames = Rle(genes$chr),
  ranges = IRanges(
    start = genes$start,
    end = genes$end,
    names = genes$gene_name
  ),
  strand = genes$strand
)


# IMPORT SALMON OUTPUT
files <- paste0("/home/artem/IMG/Projects/Elys/RNA-seq_10.07.19/",
                "larva_testes_RNA-seq_Elys_Lam_KD/salmon_out/sample",
                1:9)
files <- file.path(files, "quant.genes.sf")
names(files) <- paste0(rep(c("WT", "ElysKD", "LamKD"), each = 3), 1:3)
txi.lartest <- tximport(files, "salmon", txOut = T, importer=read.delim,
                        countsFromAbundance = "scaledTPM")

head(txi.lartest$counts)

# DESEQ2
samples <- data.frame(row.names = names(files),
                      conditions = rep(c("wt", "elysKD", "lamKD"), each = 3),
                      orig = "evrogen")

ddsTxi <- DESeqDataSetFromTximport(txi.lartest, colData = samples, design =~ conditions)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]
ddsTxi$conditions <- relevel(ddsTxi$conditions, ref = "wt")

dds <- DESeq(ddsTxi)
res.elys <- results(dds, contrast = c("conditions", "elysKD", "wt"))
res.lam <- results(dds, contrast = c("conditions", "lamKD", "wt"))
resultsNames(dds)

resLFC <- lfcShrink(dds, coef = "conditions_lamKD_vs_wt", type = "apeglm",
                    lfcThreshold = 1)
pdf("./plots/MA.plot.pdf")
DESeq2::plotMA(resLFC, alpha = 0.05, ylim=c(-6, 6))
dev.off()

summary(res.elys, alpha = 0.05)

res.df.elys <- as.data.frame(res.elys) %>% mutate(id = rownames(res.elys))
res.df.elys <- merge(genes, res.df.elys, by = "id", all.y = T)
cat((res.df.elys %>% filter(padj<0.05, log2FoldChange < -1))$id,
    file = "elys.kd.down.fbgn.txt", sep = '\n')
cat((res.df.elys %>% filter(padj<0.05, log2FoldChange > 1))$id,
    file = "elys.kd.up.fbgn.txt", sep = '\n')

res.df.lam <- as.data.frame(res.lam) %>% mutate(id = rownames(res.lam))
res.df.lam <- merge(genes, res.df.lam, by = "id", all.y = T)
cat((res.df.lam %>% filter(padj<0.05, log2FoldChange < -1))$id,
    file = "lamins.kd.down.fbgn.txt", sep = '\n')
cat((res.df.lam %>% filter(padj<0.05, log2FoldChange > 1))$id,
    file = "lamins.kd.up.fbgn.txt", sep = '\n')


save(txi.lartest, file = "RData/tximport.RData")
save(genes, res.df.elys, res.df.lam, file = "RData/lartest.DESeq.RData")
