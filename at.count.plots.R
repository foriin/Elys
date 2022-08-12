library(reshape2)
library(stringr)
library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(seqinr)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(Biostrings)
library(ggplot2)


load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/emb.hmm3.RData", verb = T)

elys.x.nup.npc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.npc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.npc <- elys.x.nup.npc[width(elys.x.nup.npc) >= 100]

elys.x.nup.nuc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.nuc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.nuc <- elys.x.nup.nuc[width(elys.x.nup.nuc) >= 100]

count_at <- function(dom){
  a <- resize(dom, 6500, fix = 'center')
  if (!all(grepl("^chr", seqlevels(a)))) seqlevels(a) <-
      paste0("chr", seqlevels(a))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, a)
  b <- sapply(as.character(seqs), function(x){
    cutt <- sapply(seq(from=1, to=nchar(x), by=100),
                   function(i) substr(x, i, i+99))
    cutt.rev <- sapply(seq(from=1, to=nchar(x), by=500),
                       function(i) substr(revcom(x), i, i+499))
    str_count(cutt, "AT") + str_count(cutt.rev, "AT")
  }, USE.NAMES = F)
  
  c <- apply(b, 1, median)
  return(c)
}

at.1 <- count_at(elys.x.nup.npc)
at.2 <- count_at(elys.x.nup.nuc)
rranges <- elys.emb.pr[start(elys.emb.pr) > 6500] %>% delete.het.gr()
at.3 <- count_at(rranges[sample(1:length(rranges), 4000)])

data <- data.frame(index = seq(-3, 3, 0.5),
                   NPC = at.1,
                   NUC = at.2,
                   random = at.3) %>% melt(id.vars = "index") %>% 
  setNames(c("kb", "Elys x", "value"))

p1 <- ggplot(data,
       aes(x = kb, y = value))+
  # geom_line()+
  # geom_line(aes(y = ymax), col = "red")+
  # geom_line(aes(y = ym), col = "blue")+
  stat_smooth(aes(col = `Elys x`),geom = "line",
              method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.2)+
  ylab("AT count per bin")+
  theme_bw()+
  ggtitle("AT-enrichment around centers of Elys domains")

pdf("plots/at_count.elys.x.nup.npc.nuc.pdf", width = 8, height = 5)
p1
dev.off()



count_at.2 <- function(dom){
  a <- resize(dom, 6020, fix = 'center')
  if (!all(grepl("^chr", seqlevels(a)))) seqlevels(a) <-
      paste0("chr", seqlevels(a))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, a)
  b <- sapply(as.character(seqs), function(x){
    cutt <- sapply(seq(from=1, to=nchar(x), by=20),
                   function(i) substr(x, i, i+19))
    1 - sapply(cutt, function(x){GC(s2c(x))}, USE.NAMES = F)
  }, USE.NAMES = F)
  
  c <- apply(b, 1, mean, na.rm = T)
  return(c)
}

at.1 <- count_at.2(elys.x.nup.npc)
at.2 <- count_at.2(elys.x.nup.nuc)
rranges <- elys.emb.pr[start(elys.emb.pr) > 6100] %>% delete.het.gr()
at.3 <- count_at.2(rranges[sample(1:length(rranges), 4000)])

data <- data.frame(index = seq(-3, 3, 0.02),
                   NPC = at.1,
                   NUC = at.2,
                   random = at.3) %>% melt(id.vars = "index") %>% 
  setNames(c("kb", "Elys x", "value"))

p2 <- ggplot(data,
             aes(x = kb, y = value))+
  # geom_line()+
  # geom_line(aes(y = ymax), col = "red")+
  # geom_line(aes(y = ym), col = "blue")+
  stat_smooth(aes(col = `Elys x`),geom = "line",
              method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.2)+
  ylab("AT-percent per 20 bp bin")+
  theme_bw()+
  ggtitle("AT-enrichment around centers of Elys domains")

pdf("plots/at_content.elys.x.nup.npc.nuc.20bp.pdf", width = 8, height = 5)
p2
dev.off()


count_at.3 <- function(dom){
  a <- resize(dom, 6100, fix = 'center')
  if (!all(grepl("^chr", seqlevels(a)))) seqlevels(a) <-
      paste0("chr", seqlevels(a))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, a)
  b <- sapply(as.character(seqs), function(x){
    cutt <- sapply(seq(from=1, to=nchar(x), by=1),
                   function(i) substr(x, i, i+99))
    1 - sapply(cutt, function(x){GC(s2c(x))}, USE.NAMES = F)
  }, USE.NAMES = F)
  
  c <- apply(b, 1, mean, na.rm = T)
  return(c)
}

at.1 <- count_at.3(elys.x.nup.npc)
at.2 <- count_at.3(elys.x.nup.nuc)
at.3 <- count_at.3(rranges[sample(1:length(rranges), 4000)])




