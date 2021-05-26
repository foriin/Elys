library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)



# Import data

load("RData/h3k27ac.RData", verbose = T)
load("RData/s2.elys.chip.RData", verbose = T)
load("RData/nup98.damid.hmm3.RData", verbose = T)
load("RData/emb.hmm.RData", verbose = T)
load("RData/emb.elys.hmm3.strict.RData", verbose = T)

seqlevels(nup98.npc.hmm3) <- sub("chr", "", seqlevels(nup98.npc.hmm3))
seqlevels(h3k27ac.hmm3) <- sub("chr", "", seqlevels(h3k27ac.hmm3))
seqlevels(nup98.nuc.hmm3) <- sub("chr", "", seqlevels(nup98.nuc.hmm3))

perm.test.jac <- function (gr1, gr2, ncores = 14) 
{
  contr.val <- bedTools.shuffle.jac(gr1, gr2, shuf = F)
  ptest <- mclapply(1:10000, function(somebodyoncetoldmetheworldsgonnarollme) {
    bedTools.shuffle.jac(gr1, gr2, shuf = T) > 
      contr.val
  }, mc.cores = ncores)
  return(sum(do.call(c, ptest))/10000)
}


nup.npc.x.h3k27ac <- sum(width(GenomicRanges::intersect(nup98.npc.hmm3, h3k27ac.hmm3, ignore.strand = T)))
nup.nuc.x.h3k27ac <- sum(width(GenomicRanges::intersect(nup98.nuc.hmm3, h3k27ac.hmm3, ignore.strand = T)))

nup.npc.x.h3k27ac.pv <- perm.test.length(nup98.npc.hmm3,
                                         h3k27ac.hmm3,
                                         contr.val = nup.npc.x.h3k27ac)
nup.npc.x.h3k27ac.pv.jac <- perm.test.jac(nup98.npc.hmm3,
                                         h3k27ac.hmm3)


nup.nuc.x.h3k27ac.pv <- perm.test.length(nup98.nuc.hmm3,
                                         h3k27ac.hmm3,
                                         contr.val = nup.nuc.x.h3k27ac)
nup.nuc.x.h3k27ac.pv.jac <- perm.test.jac(nup98.nuc.hmm3,
                                         h3k27ac.hmm3)

nup98.x.h3k27ac.df <- data.frame(NPC = c(nup.npc.x.h3k27ac,
                                         nup.npc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                         nup.npc.x.h3k27ac/sum(width(nup98.npc.hmm3))),
                                 NUC = c(nup.nuc.x.h3k27ac,
                                         nup.nuc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                         nup.nuc.x.h3k27ac/sum(width(nup98.nuc.hmm3))),
                                 row.names = c("length", "ratio from H3K27ac total",
                                               "ratio from Nup domains total"))
write.xlsx(nup98.x.h3k27ac.df, "tables/nup98.x.h3k27ac.stats.xlsx", row.names = T)

nup98.x.h3k27ac.df.2 <- tibble(NUP = c("NPC", "NUC"),
                                   Interact = c("H3K27Ac", "H3K27Ac"),
                                   ix = c(nup.npc.x.h3k27ac,
                                          nup.nuc.x.h3k27ac),
                                   `%NUP` = c(nup.npc.x.h3k27ac/sum(width(nup98.npc.hmm3)),
                                              nup.nuc.x.h3k27ac/sum(width(nup98.nuc.hmm3))),
                                   `%Interact` = c(nup.npc.x.h3k27ac/sum(width(h3k27ac.hmm3)),
                                                   nup.nuc.x.h3k27ac/sum(width(h3k27ac.hmm3))),
                                   p.v.width = c(nup.npc.x.h3k27ac.pv,
                                                 nup.nuc.x.h3k27ac.pv),
                                   p.v.jac = c(nup.npc.x.h3k27ac.pv.jac,
                                               nup.nuc.x.h3k27ac.pv.jac)
)

# Check intersection of S2 Elys ChIP peaks with Nup98 hmm3 domains
elys.x.nup.npc <- sum(width(GenomicRanges::intersect(s2.chip.bed,
                                                     nup98.npc.hmm3,
                                                     ignore.strand = T)))
elys.x.nup.nuc <- sum(width(GenomicRanges::intersect(s2.chip.bed,
                                                     nup98.nuc.hmm3,
                                                     ignore.strand = T)))
elys.x.nup.npc.pv <- perm.test.length(s2.chip.bed,
                                      nup98.npc.hmm3,
                                      contr.val = elys.x.nup.npc)
elys.x.nup.nuc.pv <- perm.test.length(s2.chip.bed,
                                      nup98.nuc.hmm3,
                                      contr.val = elys.x.nup.nuc)

# Check intersection of our embryos Elys hmm2 domains with Nup98 hmm3 domains

e.elys.x.nup.npc <- sum(width(GenomicRanges::intersect(elys.emb.hmm2,
                                                       nup98.npc.hmm3,
                                                       ignore.strand = T)))
e.elys.x.nup.nuc <- sum(width(GenomicRanges::intersect(elys.emb.hmm2,
                                                       nup98.nuc.hmm3,
                                                       ignore.strand = T)))
e.elys.x.nup.npc.pv <- perm.test.length(elys.emb.hmm2,
                                        nup98.npc.hmm3,
                                        contr.val = e.elys.x.nup.npc)
e.elys.x.nup.npc.pv.jac <- perm.test.jac(elys.emb.hmm2,
                                        nup98.npc.hmm3)
e.elys.x.nup.nuc.pv <- perm.test.length(elys.emb.hmm2,
                                        nup98.nuc.hmm3,
                                        contr.val = e.elys.x.nup.nuc)
e.elys.x.nup.nuc.pv.jac <- perm.test.jac(elys.emb.hmm2,
                                        nup98.nuc.hmm3)

# Check intersection of our embryos Elys hmm3 domains with Nup98 hmm3 domains

e.elys.x.nup.npc.2 <- sum(width(GenomicRanges::intersect(elys.emb.hmm.x,
                                                       nup98.npc.hmm3,
                                                       ignore.strand = T)))
e.elys.x.nup.nuc.2 <- sum(width(GenomicRanges::intersect(elys.emb.hmm.x,
                                                       nup98.nuc.hmm3,
                                                       ignore.strand = T)))
e.elys.x.nup.npc.2.pv <- perm.test.length(elys.emb.hmm.x,
                                        nup98.npc.hmm3,
                                        contr.val = e.elys.x.nup.npc.2)
e.elys.x.nup.npc.2.pv.jac <- perm.test.jac(elys.emb.hmm.x,
                                          nup98.npc.hmm3)
e.elys.x.nup.nuc.2.pv <- perm.test.length(elys.emb.hmm.x,
                                        nup98.nuc.hmm3,
                                        contr.val = e.elys.x.nup.nuc.2)
e.elys.x.nup.nuc.2.pv.jac <- perm.test.jac(elys.emb.hmm.x,
                                          nup98.nuc.hmm3)
elys.x.nup.df <- tibble(HMM3 = rep(c("merged", "by each rep"),each = 2),
                        NUP = rep(c("NPC", "NUC"), 2),
                        ix = c(e.elys.x.nup.npc,
                               e.elys.x.nup.nuc,
                               e.elys.x.nup.npc.2, 
                               e.elys.x.nup.nuc.2),
                        `%NUP` = c(e.elys.x.nup.npc/sum(width(nup98.npc.hmm3)),
                                   e.elys.x.nup.nuc/sum(width(nup98.nuc.hmm3)),
                                   e.elys.x.nup.npc.2/sum(width(nup98.npc.hmm3)),
                                   e.elys.x.nup.nuc.2/sum(width(nup98.nuc.hmm3))),
                        `%Elys` = c(e.elys.x.nup.npc/sum(width(elys.emb.hmm2)),
                                    e.elys.x.nup.nuc/sum(width(elys.emb.hmm2)),
                                    e.elys.x.nup.npc.2/sum(width(elys.emb.hmm.x)),
                                    e.elys.x.nup.nuc.2/sum(width(elys.emb.hmm.x))),
                        p.v.width = c(e.elys.x.nup.npc.pv,
                                      e.elys.x.nup.nuc.pv,
                                      e.elys.x.nup.npc.2.pv,
                                      e.elys.x.nup.nuc.2.pv),
                        p.v.jac = c(e.elys.x.nup.npc.pv.jac,
                                    e.elys.x.nup.nuc.pv.jac,
                                    e.elys.x.nup.npc.2.pv.jac,
                                    e.elys.x.nup.nuc.2.pv.jac)
)
save(elys.x.nup.df,
     file = "RData/emb.elys.x.nup.df.RData")
