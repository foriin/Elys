library(tibble)
library(dplyr)
library(GenomicRanges)
library(dm3)
library(rtracklayer)
library(openxlsx)

load("RData/nup98.damid.hmm3.nogap.RData", verb = T)
load("RData/s2.nups.capelson.ChIP.RData", verb = T)

# Intersection of peaks of Nup93 in S2 with hmm3 domains of Nup98 NPC* in Kc167
con.val.1 <- sum(width(GenomicRanges::intersect(s2.nup93.chip,
                                      nup98.npc.hmm3.uq.ng,
                                      ignore.strand = T)))
p.v.1 <- perm.test.length(s2.nup93.chip,
                          nup98.npc.hmm3.uq.ng,
                          contr.val = con.val.1)
# ... Nup93 peaks in S2 with Nup98 NUC* hmm3 domains
con.val.2 <- sum(width(GenomicRanges::intersect(s2.nup93.chip,
                                                nup98.nuc.hmm3.uq.ng,
                                                ignore.strand = T)))
p.v.2 <- perm.test.length(s2.nup93.chip,
                          nup98.nuc.hmm3.uq.ng,
                          contr.val = con.val.2,
                          ncores = 12)
# ... Nup98 peaks in S2 with hmm3 domains of Nup98 NPC* in Kc167
con.val.3 <- sum(width(GenomicRanges::intersect(s2.nup98.chip,
                                                nup98.npc.hmm3.uq.ng,
                                                ignore.strand = T)))
p.v.3 <- perm.test.length(s2.nup98.chip,
                          nup98.npc.hmm3.uq.ng,
                          contr.val = con.val.1, ncores = 12)
# ... Nup98 peaks in S2 with Nup98 NUC* hmm3 domains
con.val.4 <- sum(width(GenomicRanges::intersect(s2.nup98.chip,
                                                nup98.nuc.hmm3.uq.ng,
                                                ignore.strand = T)))
p.v.4 <- perm.test.length(s2.nup98.chip,
                          nup98.nuc.hmm3.uq.ng,
                          contr.val = con.val.2,
                          ncores = 12)

# Sum everything up in a table

p.v.t <- tibble(`Nup93_S2` = c(p.v.1, p.v.2),
                `Nup98_S2` = c(p.v.3, p.v.4),
                row = c("Nup98_NPC*_Kc167",
                        "Nup98_NUC*_Kc167")) %>% 
  column_to_rownames("row")

write.xlsx(p.v.t, row.names = T, "tables/nups.kc.x.nups.s2.pv.xlsx")
