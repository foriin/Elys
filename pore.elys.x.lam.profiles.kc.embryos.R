library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)

load("RData/elys.npc.x.chrom.colours.RData", verbose = T)
load("../Lam_sperm/RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
seqlevels(e.elys.npc.inactive) <- sub("chr", "", seqlevels(e.elys.npc.inactive))


# Let's make profiles of Kc167 Lamin DamID around centers of embryonic Elys
# domains that overlap with Kc167 NPC domains, LADs and inactive chromatin
# colours

# 1. Embryonic Elys domains that overlap w/ Kc NPC and LADs

domains.1 <- subsetByOverlaps(e.elys.npc.inactive[width(e.elys.npc.inactive) > 100],
                              kc.lam.hmm,
                              maxgap = 2001,
                              ignore.strand = T)
summary(width(domains.1))

domains.1.c <- resize(domains.1, 12000,fix = "center")
