library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(ggplot2)
library(gridExtra)

# Function to get damid signals around TSSs
tss.surr <- function(tss.gr, signal.gr, sig.col, vic = 5){
  ovs <- findOverlaps(tss.gr, signal.gr, ignore.strand = T)
  ovs <- ovs[ovs@to > vic]
  # print(mcols(signal.gr)[[sig.col]][1:10])
  # sapply(1:length(ovs), function(x){
  #   if(strand(tss.gr[ovs@from[x]]) %>% as.vector() == "-"){
  #     rev(mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)])
  #   } else{
  #     mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
  #   }
  # })
  
  # Don't take the strandness into account
  
  sapply(1:length(ovs), function(x){
    mcols(signal.gr)[[sig.col]][(ovs@to[x] - vic):(ovs@to[x] + vic)]
    
  })
}

rowMeans.2 <- function(df, na.rm = T){
  apply(df, 1, function(x){
    x <- x[x != 0]
    mean(x, na.rm = T)
  })
}

mean.conf <- function(vec){
  vec <- vec[vec != 0 & !is.na(vec)]
  error <- qt(0.975,df=length(vec)-1)*sd(vec)/sqrt(length(vec))
  low <- mean(vec, na.rm = T) - error
  high <- mean(vec, na.rm = T) + error
  c(low,high)
}

perm.test.ov <- function(gr.1, gr.2, shuffle = "ALL",
                         ncores = 14, N = 10000){
  cval <- length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T))
  ptest <- mclapply(1:N, function(somebodyoncetoldmetheworldsgonnarollme){
    if (shuffle %in% c("1", "ALL")) {
      gr.1 <- bedTools.shuffle.gr(gr.1)
    }
    else if (shuffle %in% c("2", "ALL")) {
      gr.2 <- bedTools.shuffle.gr(gr.2)
    }
    length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)) > cval
  }, mc.cores = ncores)
  ptest <- sum(do.call(c, ptest))
  ptest/N
}

load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("../Lam_sperm/RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/neur.enhancers.RData", verbose = T)
load("RData/starr-seq.ecd.enhancers.RData", verbose = T)

# Script used to find intersections of starr ecdysone-indused enhancers in S2 or 
# brain enhancers with NPC* or NUC* domains

# 1. Ecdysone starr enhancers
# 1.1 Intersect w/ NPC*

dom.1a <- subsetByOverlaps(nup98.npc.hmm3.uq.ng,
                                      starr.ecd,
                                      ignore.strand = T)
dom.1b <- subsetByOverlaps(nup98.nuc.hmm3.uq.ng, starr.ecd,
                           ignore.strand = T)
npc.starr.pv <- perm.test.ov(nup98.npc.hmm3.uq.ng, starr.ecd)
starr.npc.pv <- perm.test.ov(starr.ecd, nup98.npc.hmm3.uq.ng)

nuc.starr.pv <- perm.test.ov(nup98.nuc.hmm3.uq.ng, starr.ecd)
starr.nuc.pv <- perm.test.ov(starr.ecd, nup98.nuc.hmm3.uq.ng)

dom.2a <- subsetByOverlaps(nup98.npc.hmm3.uq.ng, brain.enh,
                           ignore.strand = T)
dom.2b <- subsetByOverlaps(nup98.nuc.hmm3.uq.ng, brain.enh,
                           ignore.strand = T)

npc.b.e.pv <- perm.test.ov(nup98.npc.hmm3.uq.ng, brain.enh)
nuc.b.e.pv <- perm.test.ov(nup98.nuc.hmm3.uq.ng, brain.enh)
b.e.npc.pv <- perm.test.ov(brain.enh, nup98.npc.hmm3.uq.ng)
b.e.nuc.pv <- perm.test.ov(brain.enh, nup98.nuc.hmm3.uq.ng)

