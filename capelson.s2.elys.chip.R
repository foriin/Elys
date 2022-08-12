library(dm3)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ggplot2)

load("RData/s2.elys.chip.RData", verbose = T)
elys.1 <- import.bedGraph("~/IMG/Projects/Elys/chip_capelson/covtracks/s2.elys.1.rpm.bedgraph")
elys.2 <- import.bedGraph("~/IMG/Projects/Elys/chip_capelson/covtracks/s2.elys.2.rpm.bedgraph")
input.1 <- import.bedGraph("~/IMG/Projects/Elys/chip_capelson/covtracks/s2.igg.1.rpm.bedgraph")
input.2 <- import.bedGraph("~/IMG/Projects/Elys/chip_capelson/covtracks/s2.igg.2.rpm.bedgraph")

elys.1$score <- log2(elys.1$score/input.1$score)
elys.2$score <- log2(elys.2$score/input.2$score)

export.bedGraph(elys.1[is.finite(elys.1$score)], "bed/s2.elys.chip.logpr.1.bedgraph")
export.bedGraph(elys.2[is.finite(elys.2$score)], "bed/s2.elys.chip.logpr.2.bedgraph")

elys.gr <- elys.1
elys.gr$score <- (elys.1$score + elys.2$score)/2
elys.gr <- elys.gr[is.finite(elys.gr$score)]
summary(elys.gr$score)
test <- sapply(seqlevels(elys.gr), function(x){
  summary(elys.gr$score[as.logical(seqnames(elys.gr) == x)])
})

elys.peaks.1 <- import("/home/artem/IMG/Projects/Elys/chip_capelson/peaks/elys.1_peaks.broadPeak",
                       format = "broadPeak")

elys.peaks.2 <- import("/home/artem/IMG/Projects/Elys/chip_capelson/peaks/elys.2_peaks.broadPeak",
                       format = "broadPeak")
elys.peaks <- GenomicRanges::intersect(elys.peaks.1, elys.peaks.2, ignore.strand = T)

export.bed(elys.peaks %>% renameSeqlevels(paste0("chr", seqlevels(elys.peaks))),
           con = "bed/s2.elys.peaks.bed")
export.bedGraph(elys.gr %>% renameSeqlevels(paste0("chr", seqlevels(elys.gr))),
           con = "bed/s2.elys.profile.bed")
elys.gr.2 <- subsetByOverlaps(elys.gr, elys.peaks, ignore.strand = T, maxgap = 50)

chip.log2.summ <- sapply(seqlevels(elys.gr.2), function(x){
  summary(elys.gr.2$score[as.logical(seqnames(elys.gr.2) == x)])
})

elys.df <- df.from.GRanges(elys.gr.2)
elys.df$score <- elys.gr.2$score
pdf("plots/s2.elys.log2.scores.in.peaks.pdf", width = 10, height = 4)
  boxplot(score ~ chr, data = elys.df)
  abline(h = median(elys.gr.2$score), col = "red", lty = 2)
dev.off()
