library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(openxlsx)

load("RData/nup98.damid.hmm3.strict.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)

seqlevels(nup98.npc.hmm3.uq) <- paste0("chr", seqlevels(nup98.npc.hmm3.uq))
seqlevels(nup98.nuc.hmm3.uq) <- paste0("chr", seqlevels(nup98.nuc.hmm3.uq))
seqlevels(elys.emb.hmm3) <- paste0("chr", seqlevels(elys.emb.hmm3))

colours <- read.xlsx("tables/ooo_9-state chromatin model in S2 cells.xlsx")
active.c <- colours %>% filter(Chromatin.color.type < 6) %>% makeGRangesFromDataFrame()
pc.c <- colours %>% filter(Chromatin.color.type == 6) %>% makeGRangesFromDataFrame()
silent.c <- colours %>% filter(Chromatin.color.type == 9) %>% makeGRangesFromDataFrame()

e.elys.x.nup.npc.gr <- GenomicRanges::intersect(nup98.npc.hmm3.uq,
                                                elys.emb.hmm3,
                                                ignore.strand = T)
# embryonic elys merged x nup98 npc hmm3 x chrom active states

domains.1 <- subsetByOverlaps(e.elys.x.nup.npc.gr, active.c,
                              ignore.strand = T, type = "within")

# embryonic elys merged x nup98 npc hmm3 x chrom pc bound state

domains.2 <- subsetByOverlaps(e.elys.x.nup.npc.gr, pc.c,
                              ignore.strand = T, type = "within")

# embryonic elys merged x nup98 npc hmm3 x chrom silent state

domains.3 <- subsetByOverlaps(e.elys.x.nup.npc.gr, silent.c,
                              ignore.strand = T, type = "within")

export.bed(domains.1, "bed/nup.npc.x.elys.merged.active.bed")
export.bed(domains.2, "bed/nup.npc.x.elys.merged.PC.bed")
export.bed(domains.3, "bed/nup.npc.x.elys.merged.silent.bed")
export.bed(active.c, "bed/active.chrom.colours.bed")
