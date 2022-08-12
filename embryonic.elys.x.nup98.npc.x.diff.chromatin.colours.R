library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(openxlsx)

load("RData/nup98.damid.hmm3.strict.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)

# Add "chr" to chromosome names for exporting in UCSC format
seqlevels(nup98.npc.hmm3.uq) <- paste0("chr", seqlevels(nup98.npc.hmm3.uq))
seqlevels(nup98.nuc.hmm3.uq) <- paste0("chr", seqlevels(nup98.nuc.hmm3.uq))
seqlevels(elys.emb.hmm3) <- paste0("chr", seqlevels(elys.emb.hmm3))
# Close 900 nt gaps
elys.emb.hmm3 <- reduce(elys.emb.hmm3, min.gapwidth = 901)
# Load chromatin colours data (Kharchenko et al., 2011)
colours <- read.xlsx("tables/ooo_9-state chromatin model in S2 cells.xlsx")
col.gr <- makeGRangesFromDataFrame(colours %>% 
                           mutate(Chr = sub("chr", "", Chr)) %>% 
                             rename("Chromatin.color.type" = "colour") %>% 
                           select(1:4), keep.extra.columns = T)
# Make Granges of active chromatin states (colours 1-6 in Kharchenko article)
active.c <- colours %>% filter(Chromatin.color.type < 6) %>%
  makeGRangesFromDataFrame() %>% reduce(min.gapwidth=100)
# Make GRanges of inactive chromatin states (states 6-9) in Kharchenko article
inactive.c <- colours %>% filter(Chromatin.color.type %in% c(6,7,8,9)) %>%
  makeGRangesFromDataFrame() %>% reduce()
# Intersection of embryonic hmm3 Elys w/ Nup98 NPC non-intersecting w/ Nup98 Nuc
e.elys.x.nup.npc.gr <- GenomicRanges::intersect(nup98.npc.hmm3.uq,
                                                elys.emb.hmm3,
                                                ignore.strand = T)

# embryonic elys merged x nup98 npc hmm3 x chrom active states
e.elys.npc.active <- subsetByOverlaps(e.elys.x.nup.npc.gr, active.c,
                              ignore.strand = T, type = "within") %>% 
  reduce(min.gapwidth = 901)

# embryonic elys merged x nup98 npc hmm3 x chrom inactive states
e.elys.npc.inactive <- subsetByOverlaps(e.elys.x.nup.npc.gr, inactive.c,
                              ignore.strand = T, type = "within") %>% 
  reduce(min.gapwidth = 901)

# Export to .bed
export.bed(e.elys.npc.active, "bed/nup.npc.x.elys.merged.active.bed")
export.bed(e.elys.npc.inactive, "bed/nup.npc.x.elys.merged.inactive.bed")
export.bed(active.c, "bed/active.chrom.colours.r.bed")

# remove "chr" so that you could work with it later
seqlevels(nup98.npc.hmm3.uq) <- sub("chr", "", seqlevels(nup98.npc.hmm3.uq))
seqlevels(nup98.nuc.hmm3.uq) <- sub("chr", "", seqlevels(nup98.nuc.hmm3.uq))
seqlevels(elys.emb.hmm3) <- sub("chr", "", seqlevels(elys.emb.hmm3))
# save for subsequent use
save(e.elys.npc.active, e.elys.npc.inactive,
     e.elys.x.nup.npc.gr, file = "RData/elys.npc.x.chrom.colours.RData")
save(col.gr, file = "RData/kharchenko.chromatin.colours.RData")
