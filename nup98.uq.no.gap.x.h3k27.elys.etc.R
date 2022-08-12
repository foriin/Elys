library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)



# Import data

load("RData/h3k27ac.RData", verbose = T)
load("RData/s2.elys.chip.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)

seqlevels(h3k27ac.hmm3, pruning.mode = "coarse") <- euc.chroms
seqlevels(elys.emb.hmm3.ng, pruning.mode = "coarse") <- euc.chroms

perm.test.jac <- function (gr1, gr2, ncores = 14) 
{
  contr.val <- bedTools.shuffle.jac(gr1, gr2, shuf = F)
  ptest <- mclapply(1:10000, function(somebodyoncetoldmetheworldsgonnarollme) {
    bedTools.shuffle.jac(gr1, gr2, shuf = T) > 
      contr.val
  }, mc.cores = ncores)
  return(sum(do.call(c, ptest))/10000)
}

ndom.overlap <- function(gr1, gr2){
  a <- length(subsetByOverlaps(gr1, gr2, ignore.strand = T))/length(gr1)
  b <- length(subsetByOverlaps(gr2, gr1, ignore.strand = T))/length(gr2)
  c(a,b)
}
nup.npc.x.h3k27ac <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                                        h3k27ac.hmm3,
                                                        ignore.strand = T)))
nup.nuc.x.h3k27ac <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3.uq.ng,
                                                        h3k27ac.hmm3,
                                                        ignore.strand = T)))

nup.npc.x.h3k27ac.pv <- perm.test.length(nup98.npc.hmm3.uq.ng,
                                         h3k27ac.hmm3,
                                         contr.val = nup.npc.x.h3k27ac)
# nup.npc.x.h3k27ac.pv.jac <- perm.test.jac(nup98.npc.hmm3.uq,
# h3k27ac.hmm3)


nup.nuc.x.h3k27ac.pv <- perm.test.length(nup98.nuc.hmm3.uq.ng,
                                         h3k27ac.hmm3,
                                         contr.val = nup.nuc.x.h3k27ac)
ndom.overlap(nup98.nuc.hmm3.uq.ng, h3k27ac.hmm3)
# nup.nuc.x.h3k27ac.pv.jac <- perm.test.jac(nup98.nuc.hmm3.uq,
# h3k27ac.hmm3)

nup98.x.h3k27ac.df <- data.frame(NPC = c(nup.npc.x.h3k27ac,
                                         nup.npc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                         nup.npc.x.h3k27ac/sum(width(nup98.npc.hmm3.uq.ng))),
                                 NUC = c(nup.nuc.x.h3k27ac,
                                         nup.nuc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                         nup.nuc.x.h3k27ac/sum(width(nup98.nuc.hmm3.uq.ng))),
                                 row.names = c("length", "ratio from H3K27ac total",
                                               "ratio from Nup domains total"))
# write.xlsx(nup98.x.h3k27ac.df, "tables/nup98.x.h3k27ac.stats.xlsx", row.names = T)

nup98.x.h3k27ac.df.2 <- tibble(NUP = c("NPC*", "NUC*"),
                               Interact = c("H3K27Ac", "H3K27Ac"),
                               ix = c(nup.npc.x.h3k27ac,
                                      nup.nuc.x.h3k27ac),
                               `%NUP length` = c(nup.npc.x.h3k27ac/sum(width(nup98.npc.hmm3.uq.ng)),
                                                 nup.nuc.x.h3k27ac/sum(width(nup98.nuc.hmm3.uq.ng))),
                               `%Interact length` = c(nup.npc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                                      nup.nuc.x.h3k27ac/sum(width(h3k27ac.hmm3))),
                               p.v.width = c(nup.npc.x.h3k27ac.pv,
                                             nup.nuc.x.h3k27ac.pv)
) %>% cbind(rbind(ndom.overlap(nup98.npc.hmm3.uq.ng, h3k27ac.hmm3),
                  ndom.overlap(nup98.nuc.hmm3.uq.ng, h3k27ac.hmm3)) %>% 
              as.data.frame() %>% setNames(c("%NUP domains N", "%Interact domains N")))

# Check intersection of S2 Elys ChIP peaks with Nup98 hmm3 domains
elys.x.nup.npc <- sum(width(GenomicRanges::intersect(s2.chip.bed,
                                                     nup98.npc.hmm3.uq.ng,
                                                     ignore.strand = T)))
elys.x.nup.nuc <- sum(width(GenomicRanges::intersect(s2.chip.bed,
                                                     nup98.nuc.hmm3.uq.ng,
                                                     ignore.strand = T)))
elys.x.nup.npc.pv <- perm.test.length(s2.chip.bed,
                                      nup98.npc.hmm3.uq.ng,
                                      contr.val = elys.x.nup.npc)
elys.x.nup.nuc.pv <- perm.test.length(s2.chip.bed,
                                      nup98.nuc.hmm3.uq.ng,
                                      contr.val = elys.x.nup.nuc)

# Check intersection of our embryos Elys merged hmm3 domains with Nup98 hmm3 domains

e.elys.x.nup.npc <- sum(width(GenomicRanges::intersect(elys.emb.hmm3.ng,
                                                       nup98.npc.hmm3.uq.ng,
                                                       ignore.strand = T)))
e.elys.x.nup.nuc <- sum(width(GenomicRanges::intersect(elys.emb.hmm3.ng,
                                                       nup98.nuc.hmm3.uq.ng,
                                                       ignore.strand = T)))
e.elys.x.nup.npc.pv <- perm.test.length(elys.emb.hmm3.ng,
                                        nup98.npc.hmm3.uq.ng,
                                        contr.val = e.elys.x.nup.npc)
# e.elys.x.nup.npc.pv.jac <- perm.test.jac(elys.emb.hmm3,
# nup98.npc.hmm3.uq)
e.elys.x.nup.nuc.pv <- perm.test.length(elys.emb.hmm3.ng,
                                        nup98.nuc.hmm3.uq.ng,
                                        contr.val = e.elys.x.nup.nuc)
# e.elys.x.nup.nuc.pv.jac <- perm.test.jac(elys.emb.hmm3,
# nup98.nuc.hmm3.uq)
e.elys.x.nup.npc.nd <- length(subsetByOverlaps(elys.emb.hmm3, nup98.npc.hmm3.uq,
                                               ignore.strand = T))

elys.x.nup.uq.df <- tibble(Tissue = rep(c("S2", "Embryo"), each = 2),
                           HMM3 = rep(c("NA", "merged"),each = 2),
                           NUP = rep(c("NPC*", "NUC*"), 2),
                           ix = c(elys.x.nup.npc,
                                  elys.x.nup.nuc,
                                  e.elys.x.nup.npc,
                                  e.elys.x.nup.nuc),
                           `%NUP length` = c(elys.x.nup.npc/sum(width(nup98.npc.hmm3.uq.ng)),
                                             elys.x.nup.nuc/sum(width(nup98.nuc.hmm3.uq.ng)),
                                             e.elys.x.nup.npc/sum(width(nup98.npc.hmm3.uq.ng)),
                                             e.elys.x.nup.nuc/sum(width(nup98.nuc.hmm3.uq.ng))),
                           `%Elys length` = c(elys.x.nup.npc/sum(width(s2.chip.bed)),
                                              elys.x.nup.nuc/sum(width(s2.chip.bed)),
                                              e.elys.x.nup.npc/sum(width(elys.emb.hmm3.ng)),
                                              e.elys.x.nup.nuc/sum(width(elys.emb.hmm3.ng))),
                           p.v.width = c(elys.x.nup.npc.pv,
                                         elys.x.nup.nuc.pv,
                                         e.elys.x.nup.npc.pv,
                                         e.elys.x.nup.nuc.pv)) %>% 
  cbind(rbind(
    ndom.overlap(nup98.npc.hmm3.uq.ng, s2.chip.bed),
    ndom.overlap(nup98.nuc.hmm3.uq.ng, s2.chip.bed),
    ndom.overlap(nup98.npc.hmm3.uq.ng, elys.emb.hmm3.ng),
    ndom.overlap(nup98.nuc.hmm3.uq.ng, elys.emb.hmm3.ng)
  ) %>% as.data.frame() %>% setNames(c(
    "%NUP domains N", "%Elys domains N"
  ))
  )
save(elys.x.nup.uq.df,
     file = "RData/emb.elys.x.nup.uq.ng.df.RData")
