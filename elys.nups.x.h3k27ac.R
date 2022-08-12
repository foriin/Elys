library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(ggplot2)
library(gridExtra)
library(openxlsx)



jacc.ov <- function(gr1,gr2){
  ov <- findOverlaps(gr1,gr2)
  jac = (unique(ov@from) %>% length + unique(ov@to) %>% length)/(length(gr1) + length(gr2))
  jac
}

perm.test.ov <- function(gr.1, gr.2, shuffle = "ALL",
                         ncores = 14, N = 10000,
                         genome = "/home/artem/Database/dmel/Genome/dm3.genome"){
  cval.1 <- length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T))
  cval.2 <- length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T))
  ptest <- mclapply(1:N, function(somebodyoncetoldmetheworldsgonnarollme){
    if (shuffle == "ALL") {
      gr.1 <- bedTools.shuffle.gr(gr.1, genome = genome)
      gr.2 <- bedTools.shuffle.gr(gr.2, genome = genome)
    }
    else if (shuffle %in% c("1", "2")) {
      assign(paste0("gr.", shuffle),
             bedTools.shuffle.gr(get(paste0("gr.", shuffle))),
             genome = genome)
    }
    c(length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
      length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)) > cval.1,
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)) > cval.2)
  }, mc.cores = ncores)
  ptest <- do.call(rbind, ptest)
  pv.1 <- sum(ptest[,3])/N
  pv.2 <- sum(ptest[,4])/N
  avg.rand.len.1 <- mean(ptest[,1])
  avg.rand.len.2 <- mean(ptest[,2])
  c("AxB" = cval.1,
    "AxB random" = avg.rand.len.1,
    "%A" = cval.1/length(gr.1)*100,
    "BxA" = cval.2,
    "BxA random" = avg.rand.len.2,
    "%B" = cval.2/length(gr.2)*100,
    "Pval A" = pv.1,
    "Pval B" = pv.2)
}
perm.test.jac.ov <- function(gr.1, gr.2, shuffle = "ALL",
                             ncores = 14, N = 10000,
                             genome = "/home/artem/Database/dmel/Genome/dm3.genome"){
  cval <- jacc.ov(gr.1, gr.2)
  ptest <- mclapply(1:N, function(somebodyoncetoldmetheworldsgonnarollme){
    if (shuffle == "ALL") {
      gr.1 <- bedTools.shuffle.gr(gr.1, genome = genome)
      gr.2 <- bedTools.shuffle.gr(gr.2, genome = genome)
    }
    else if (shuffle %in% c("1", "2")) {
      assign(paste0("gr.", shuffle),
             bedTools.shuffle.gr(get(paste0("gr.", shuffle))),
             genome = genome)
    }
    c(length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
      jacc.ov(gr.1, gr.2) > cval)
  }, mc.cores = ncores)
  ptest <- do.call(rbind, ptest)
  pv.1 <- sum(ptest[,3])/N
  avg.rand.len.1 <- mean(ptest[,1])
  avg.rand.len.2 <- mean(ptest[,2])
  c("AxB" = length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
    "AxB random" = avg.rand.len.1,
    "%A" = length(subsetByOverlaps(gr.1, gr.2,
                                   ignore.strand = T))/length(gr.1)*100,
    "BxA" = length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
    "BxA random" = avg.rand.len.2,
    "%B" =length(subsetByOverlaps(gr.2, gr.1,
                                  ignore.strand = T))/length(gr.2)*100,
    "Pval" = pv.1)
}

load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/h3k27ac.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/nup98.damid.hmm3.RData", verb = T)
load("RData/neur.enhancers.RData", verbose = T)
load("RData/starr-seq.ecd.enhancers.RData", verbose = T)
load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/brain.enhancers.RData", verb = T)
load("RData/s2.elys.difex.tpm.RData", verbose = T)
load("RData/kharchenko.chromatin.colours.RData", verbose = T)


# NUP98 domains that are common between NPC and nucleoplasmic

nup98.dual <- GenomicRanges::intersect(nup98.npc.hmm3,
                                       nup98.nuc.hmm3,
                                       ignore.strand = T)
nup98.dual.x.elys.emb <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                                  nup98.dual,
                                                  ignore.strand = T)
export.bed(nup98.dual.x.elys.emb %>%
             renameSeqlevels(paste0("chr", seqlevels(nup98.dual.x.elys.emb))),
           "bed/elys.emb.hmm3.x.nup98.dual.bed")


e.npc.h3k27ac.pt.jac <- perm.test.jac.ov(nup.npc.x.elys.emb.1,
                                         h3k27ac.hmm3)
e.nuc.h3k27ac.pt.jac <- perm.test.jac.ov(nup.nuc.x.elys.emb.1,
                                         h3k27ac.hmm3)
e.dual.h3k27ac.pt.jac <- perm.test.jac.ov(nup98.dual.x.elys.emb,
                                          h3k27ac.hmm3)

e.npc.h3k27ac.pt.ov <- perm.test.ov(nup.npc.x.elys.emb.1,
                                         h3k27ac.hmm3)
e.nuc.h3k27ac.pt.ov <- perm.test.ov(nup.nuc.x.elys.emb.1,
                                         h3k27ac.hmm3)
e.dual.h3k27ac.pt.ov <- perm.test.ov(nup98.dual.x.elys.emb,
                                          h3k27ac.hmm3)
tab.h3k27ac <- cbind(
  tibble("A" = c("Elys x Nup NPC", "Elys x Nup NUC", "Elys x Nup Dual"),
         "B" = rep("S2 H3K27Ac HMM3", 3)),
  rbind(e.npc.h3k27ac.pt.ov,
                     e.nuc.h3k27ac.pt.ov,
                     e.dual.h3k27ac.pt.ov)) %>% as_tibble() %>% 
  mutate("Pv jac" = c(e.npc.h3k27ac.pt.jac[7],
                      e.nuc.h3k27ac.pt.jac[7],
                      e.dual.h3k27ac.pt.jac[7]))
write.xlsx(tab.h3k27ac, "tables/elys.x.nup.h3k27ac.xlsx")

b.e.900.genes <- res.df.elyskd %>% filter(id %in% b.e.genes) 
write.xlsx(b.e.900.genes, "tables/s2.rnaseq.900.enh.genes.xlsx")




