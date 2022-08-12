library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)

# Here I'll only use npc and nucleoplasmic Nup98 domains that do not overlap
# with each other. Also, no jaccard tests
# UPD. We filled the gaps lesser than 900 bp in Nup domains mentioned above
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/s2.gaf.and.pbap.bap.chip.peaks.RData", verbose = T)


ndom.overlap <- function(gr1, gr2){
  a <- length(subsetByOverlaps(gr1, gr2, ignore.strand = T))/length(gr1)
  b <- length(subsetByOverlaps(gr2, gr1, ignore.strand = T))/length(gr2)
  c(a,b)
}


nup.npc.x.gaf <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                                    s2.gaf.conf.peaks,
                                                    ignore.strand = T)))
nup.nuc.x.gaf <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3.uq.ng,
                                                    s2.gaf.conf.peaks,
                                                    ignore.strand = T)))

nup.npc.x.gaf.pv <- perm.test.length(nup98.npc.hmm3.uq.ng,
                                      s2.gaf.conf.peaks %>% sort,
                                     contr.val = nup.npc.x.gaf)

nup.nuc.x.gaf.pv <- perm.test.length(nup98.nuc.hmm3.uq.ng,
                                      s2.gaf.conf.peaks %>% sort,
                                  contr.val = nup.nuc.x.gaf)

nup.npc.x.pbap.bap <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                                         s2.pbap.bap.peaks,
                                                         ignore.strand = T)))
nup.nuc.x.pbap.bap <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3.uq.ng,
                                                         s2.pbap.bap.peaks,
                                                         ignore.strand = T)))

nup.npc.x.pbap.bap.pv <- perm.test.length(nup98.npc.hmm3.uq.ng,
                                           s2.pbap.bap.peaks,
                                          contr.val = nup.npc.x.pbap.bap)

nup.nuc.x.pbap.bap.pv <- perm.test.length(nup98.nuc.hmm3.uq.ng,
                                           s2.pbap.bap.peaks,
                                          contr.val = nup.nuc.x.pbap.bap)
nup98.x.df.uq <- tibble(NUP = rep(c("NPC*", "NUC*"),2),
                               Interact = rep(c("GAF", "pBAP/BAP"), each = 2),
                               ix = c(nup.npc.x.gaf,
                                      nup.nuc.x.gaf,
                                      nup.npc.x.pbap.bap,
                                      nup.nuc.x.pbap.bap),
                               `%NUP length` = c(nup.npc.x.gaf/sum(width(nup98.npc.hmm3.uq.ng)),
                                          nup.nuc.x.gaf/sum(width(nup98.nuc.hmm3.uq.ng)),
                                          nup.npc.x.pbap.bap/sum(width(nup98.npc.hmm3.uq.ng)),
                                          nup.nuc.x.pbap.bap/sum(width(nup98.nuc.hmm3.uq.ng))),
                               `%Interact length` = c(nup.npc.x.gaf/sum(width(s2.gaf.conf.peaks)),
                                               nup.nuc.x.gaf/sum(width(s2.gaf.conf.peaks)),
                                               nup.npc.x.pbap.bap/sum(width(s2.pbap.bap.peaks)),
                                               nup.nuc.x.pbap.bap/sum(width(s2.pbap.bap.peaks))),
                               p.v.width = c(nup.npc.x.gaf.pv,
                                             nup.nuc.x.gaf.pv,
                                             nup.npc.x.pbap.bap.pv,
                                             nup.nuc.x.pbap.bap.pv)
) %>% cbind(rbind(ndom.overlap(nup98.npc.hmm3.uq.ng, s2.gaf.conf.peaks),
                  ndom.overlap(nup98.nuc.hmm3.uq.ng, s2.gaf.conf.peaks),
                  ndom.overlap(nup98.npc.hmm3.uq.ng, s2.pbap.bap.peaks),
                  ndom.overlap(nup98.nuc.hmm3.uq.ng, s2.pbap.bap.peaks)
                  ) %>% as.data.frame() %>% 
              setNames(c("%NUP domains N", "%Interact domains N")))

nup98.x.df.uq <- rbind(nup98.x.h3k27ac.df.2, nup98.x.df.uq)
save(nup98.x.df.uq,
     file = "RData/nup98.uq.ng.dom.ix.table.RData")
