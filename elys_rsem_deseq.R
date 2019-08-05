library(GenomicRanges)
library(dplyr)
library(edgeR)
library(DESeq2)
library(tximport)
library(data.table)


rm(list = ls())

setwd("~/IMG/Projects/Elys/kd_rnaseq/")

emb.cnt <- fread("Emb_control_rsem.genes.results")
emb.kde <- fread("Emb_kdElys_rsem.genes.results")

s2.cnt <- fread("S2_control_rsem.genes.results")
s2.kde <- fread("S2_kdElys_rsem.genes.results")

rsem.out <- tximport(dir(pattern = "S2.*genes"),
                     type = "rsem")

samples <- data.frame(row.names = c("S2_cnt", "S2_kd"), conditions = c("wt", "kd"))
rsem.out$length[which(rsem.out$length == 0)] <- 1

ddsTxi <- DESeqDataSetFromTximport(rsem.out, colData = samples, design =~ conditions)
counts(ddsTxi)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]

dds <- DESeq(ddsTxi)
res <- results(dds)
resultsNames(dds)

resOrdered <- res[order(res$padj),]
summary(res)

res.df <- as.data.frame(resOrdered) %>% mutate(fbgn = sub("_.*", "", row.names(resOrdered)),
                                               gene_name = sub("[^_]+_([^:]+):?.*", "\\1", row.names(resOrdered))) %>% filter(!is.na(padj))

##############################################################

rsem.out <- tximport(dir(pattern = "Emb.*genes"),
                     type = "rsem")

samples <- data.frame(row.names = c("emb_cnt", "emb_kd"), conditions = c("wt", "kd"))
rsem.out$length[which(rsem.out$length == 0)] <- 1

ddsTxi <- DESeqDataSetFromTximport(rsem.out, colData = samples, design =~ conditions)
counts(ddsTxi)

keep <- rowSums(counts(ddsTxi)) >= 10

ddsTxi <- ddsTxi[keep,]

dds <- DESeq(ddsTxi)
res <- results(dds)
resultsNames(dds)

resOrdered <- res[order(res$padj),]
summary(res)

res.df <- as.data.frame(resOrdered) %>% mutate(fbgn = sub("_.*", "", row.names(resOrdered)),
                                               gene_name = sub("[^_]+_([^:]+):?.*", "\\1", row.names(resOrdered))) %>% filter(!is.na(padj))


