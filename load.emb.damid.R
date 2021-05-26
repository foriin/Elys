library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)


lam.emb.pr <- read.table("~/IMG/Projects/Elys/26.03.18_elys_lam_embryos/500nt/Bedgraph/LAM.EMB.wt.sum.norm.bedgraph",
                         skip = 1) %>% mutate(V1 = sub("chr", "", V1),
                                              V2 = V2 + 1) %>% 
  setNames(c("chr", "start", "end", "score")) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

elys.emb.pr <- read.table("~/IMG/Projects/Elys/26.03.18_elys_lam_embryos/500nt/Bedgraph/ELYS.EMB.wt.sum.norm.bedgraph",
                         skip = 1) %>% mutate(V1 = sub("chr", "", V1),
                                              V2 = V2 + 1) %>% 
  setNames(c("chr", "start", "end", "score")) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

lam.emb.hmm <- import.bed("~/IMG/Projects/Elys/26.03.18_elys_lam_embryos/500nt/BioHMM/LAM.EMB.500nt.domains.bed") %>% 
  GRanges()
elys.emb.hmm <- import.bed("~/IMG/Projects/Elys/26.03.18_elys_lam_embryos/500nt/BioHMM/ELYS.EMB.500nt.domains.bed") %>% 
  GRanges()
summary(width(lam.emb.hmm))
summary(width(elys.emb.hmm))
summary(lam.emb.pr$score)
summary(elys.emb.pr$score)
