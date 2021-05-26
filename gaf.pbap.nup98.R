library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)

chroms <- c("2L", "2R", "3L", "3R", "X")
perm.test.jac <- function (gr1, gr2, ncores = 14) 
{
  contr.val <- bedTools.shuffle.jac(gr1, gr2, shuf = F)
  ptest <- mclapply(1:10000, function(somebodyoncetoldmetheworldsgonnarollme) {
      bedTools.shuffle.jac(gr1, gr2, shuf = T) > 
      contr.val
  }, mc.cores = ncores)
  return(sum(do.call(c, ptest))/10000)
}

load("RData/nup98.damid.hmm3.RData", verbose = T)

s2.gaf.conf.peaks <- import.bed("bed/s2.gaf.conf.peaks.bed")
s2.gaf.conf.peaks <- s2.gaf.conf.peaks[seqnames(s2.gaf.conf.peaks) %in% chroms]
seqlevels(s2.gaf.conf.peaks) <- seqlevelsInUse(s2.gaf.conf.peaks)
seqlevels(s2.gaf.conf.peaks) <- chroms
# seqlevels(s2.gaf.conf.peaks) <- sub("chr", "", seqlevels(s2.gaf.conf.peaks))

s2.pbap.bap.peaks <- import.bed("bed/pbap.bap.bed")
seqlevels(s2.pbap.bap.peaks) <- sub("chr", "", seqlevels(s2.pbap.bap.peaks))

nup.npc.x.gaf <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3, s2.gaf.conf.peaks, ignore.strand = T)))
nup.nuc.x.gaf <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3, s2.gaf.conf.peaks, ignore.strand = T)))

nup.npc.x.gaf.pv <- perm.test.length(nup98.npc.hmm3,
                                      s2.gaf.conf.peaks %>% sort,
                                     contr.val = nup.npc.x.gaf)
nup.npc.x.gaf.pv.jac <- perm.test.jac(nup98.npc.hmm3,
                                         s2.gaf.conf.peaks %>% sort)

nup.nuc.x.gaf.pv <- perm.test.length(nup98.nuc.hmm3,
                                      s2.gaf.conf.peaks %>% sort,
                                  contr.val = nup.nuc.x.gaf)
nup.nuc.x.gaf.pv.jac <- perm.test.jac(nup98.nuc.hmm3,
                                         s2.gaf.conf.peaks %>% sort)

nup.npc.x.pbap.bap <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3, s2.pbap.bap.peaks, ignore.strand = T)))
nup.nuc.x.pbap.bap <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3, s2.pbap.bap.peaks, ignore.strand = T)))

nup.npc.x.pbap.bap.pv.jac <- perm.test.jac(nup98.npc.hmm3,
                                     s2.pbap.bap.peaks)
nup.npc.x.pbap.bap.pv <- perm.test.length(nup98.npc.hmm3,
                                           s2.pbap.bap.peaks,
                                          contr.val = nup.npc.x.pbap.bap)

nup.nuc.x.pbap.bap.pv.jac <- perm.test.jac(nup98.nuc.hmm3,
                                     s2.pbap.bap.peaks)
nup.nuc.x.pbap.bap.pv <- perm.test.length(nup98.nuc.hmm3,
                                           s2.pbap.bap.peaks,
                                          contr.val = nup.nuc.x.pbap.bap)
nup98.x.df.2 <- tibble(NUP = rep(c("NPC", "NUC"),2),
                               Interact = rep(c("GAF", "pBAP/BAP"), each = 2),
                               ix = c(nup.npc.x.gaf,
                                      nup.nuc.x.gaf,
                                      nup.npc.x.pbap.bap,
                                      nup.nuc.x.pbap.bap),
                               `%NUP` = c(nup.npc.x.gaf/sum(width(nup98.npc.hmm3)),
                                          nup.nuc.x.gaf/sum(width(nup98.nuc.hmm3)),
                                          nup.npc.x.pbap.bap/sum(width(nup98.npc.hmm3)),
                                          nup.nuc.x.pbap.bap/sum(width(nup98.nuc.hmm3))),
                               `%Interact` = c(nup.npc.x.gaf/sum(width(s2.gaf.conf.peaks)),
                                               nup.nuc.x.gaf/sum(width(s2.gaf.conf.peaks)),
                                               nup.npc.x.pbap.bap/sum(width(s2.pbap.bap.peaks)),
                                               nup.nuc.x.pbap.bap/sum(width(s2.pbap.bap.peaks))),
                               p.v.width = c(nup.npc.x.gaf.pv,
                                             nup.nuc.x.gaf.pv,
                                             nup.npc.x.pbap.bap.pv,
                                             nup.nuc.x.pbap.bap.pv),
                               p.v.jac = c(nup.npc.x.gaf.pv.jac,
                                           nup.nuc.x.gaf.pv.jac,
                                           nup.npc.x.pbap.bap.pv.jac,
                                           nup.nuc.x.pbap.bap.pv.jac)
)

nup.98.x.df <- rbind(nup98.x.h3k27ac.df.2, nup98.x.df.2)
save(nup.98.x.df,
     file = "RData/nup98.dom.ix.table.RData")
