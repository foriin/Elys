library(dplyr)
library(GenomicRanges)
library(snapCGH)
.libPaths(c(.libPaths(), "/home/artem/R/x86_64-pc-linux-gnu-library/3.5/"))

rm(list = ls())

# 
# pc.chip <- read.table("/home/artem/IMG/Projects/LAM.SpG.SpC/Testes_Pc_ChIP/GPL15057-6120 Запрос.txt")
# pc.chip <- pc.chip %>% select(4,5,6,9)
# unique(pc.chip$V4)
load("Pc_chip.RData")

source("/home/artem/R/projects/DamIDbinscount/DamID_HMM.functions.R")


pc.chip <- cbind(paste0(pc.chip$V4, 1:nrow(pc.chip)), pc.chip) %>%
  setNames(c("ID", "chr", "start", "end", "DamID.value"))

pc.chip.sp <- split(pc.chip, pc.chip$chr)
lastcol <- 3
list(pc.chip.sp)[[1]][[1]] %>% head
runBioHMM.2
fit.model.2
hmm.2.wrapper(DATAs)

DATA <- sooka %>% setNames(c("chr", "start", "end", "DamID.value"))

lastcol <- 3

DATAs <- lapply(DATA[(lastcol + 1):ncol(DATA)], function(x){
  names(x) = lapply(chroms, function(y) filter(cbind(DATA[1:lastcol],
                                                     "DamID.value" = x), chr == y,
                                               !is.na(DamID.value)))
})

DATAs[[1]]
str(pc.chip)
chroms <- paste0("chr", chroms)

zoop <- hmm.2.wrapper(DATAs)
pc.chip.hmm.3 <- hmm.3.wrapper(DATAs)
pc.chip.hmm.3.gr <- makeGRangesFromDataFrame(pc.chip.hmm.3[[1]] %>%
                                               filter(domain == 1)) %>% 
  GenomicRanges::reduce()

export.bed(pc.chip.hmm.3.gr, "pc.testes.hmm.3.bed")


oof <- runBioHMM.2(DATAs[[1]][[1]]$DamID.value, DATAs[[1]][[1]][, -5])
library(rtracklayer)

write.table(pc.chip[, -1], "fuckthisshit.bedgraph", quote = F, sep = "\t", row.names = F,
            col.names = F)

sooka <- fread("oof.bed") %>% filter(V4 != ".")
sooka$V4 <- as.numeric(sooka$V4)

zeep <- makeGRangesFromDataFrame(zoop[[1]] %>% filter(domain == 1))
zeep <- GenomicRanges::reduce(zeep)

export.bed(zeep, "pc.testes.hmm2.bed")
