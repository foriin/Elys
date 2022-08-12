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

# Metagene plots of AT-percent for genes that overlap w/ Elys x Nup98 NPC,
# Elys x Nup98 nucleoplasmic or Elys domains that don't overlap with Nup98
# Here I used 25 bp windows for counting AT percent


load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/emb.hmm3.RData", verb = T)
load("RData/nups.x.elys.emb.RData", verbose = T)


bins <- GRanges(elys.emb.pr)
seqlevels(bins) <- paste0("chr", seqlevels(bins))
bins.seq <- as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm3,bins))
score <- mclapply(bins.seq, function(string){
  1 - GC(s2c(as.character(string)))
}, mc.cores = 12)
bins$score <- unlist(score)
seqlevels(bins) <- sub("chr", "", seqlevels(bins))

bw.pre <- bins[!is.na(bins$score)]
seqlengths(bw.pre) <- chr.lengths[match(names(seqlengths(bw.pre)), names(chr.lengths))]
export.bw(bw.pre, "metagene/dm3.300bp.at.perc.bigWig")
# import deeptools data
at.npc <- read.table("metagene/metaplot.elys.npc.at.count.300bp.tab", skip = 1)[, -1] %>% t()
at.nuc <- read.table("metagene/metaplot.elys.nuc.at.count.300bp.tab", skip = 1)[, -1] %>% t()
at.oth <- read.table("metagene/metaplot.elys.other.at.count.300bp.tab", skip = 1)[, -1] %>% t()
at.pr <- as.data.frame(cbind(at.npc, at.nuc, at.oth))
at.pr <- at.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

at.m <- melt(at.pr, id.vars = "kb")

p.at <- ggplot(at.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("blue",
                                "green",
                                "purple"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("% AT")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/o Nups")

pdf("plots/metagene.AT.percent.300bp.pdf")
p.at
dev.off()

load("RData/elys.nup.npc.metagene.RData", verbose = T)
elys.m <- melt(elys.pr, id.vars = "kb")
elys.m <- elys.m %>% setNames(c("kb", "variable", "elys"))
at.m <- at.m %>% setNames(c("kb", "variable", "AT"))
elys.at <- merge(elys.m, at.m, by = c("kb", "variable"))

at.labs <-  seq(0.4, 0.9, by = 0.05)
at.scaled <- (at.labs * 3) - 1.8

p1 <- ggplot(elys.at %>% 
         filter(variable == "npc"), aes(x = kb))+
  geom_line(aes(y = elys),col = "blue")+
  geom_line(aes(y = AT * 3 - 1.8), col = "red")+
  scale_y_continuous(
    name = "log2(Dam-Elys/Dam)",
    limits = c(-0.5, 0.9),
    sec.axis = sec_axis(trans = (~.),
                        breaks = at.scaled,
                        labels = at.labs,
                        name = "% AT")
  )+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "blue", size=13),
        axis.title.y.right = element_text(color = "red", size=13))+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  # ylim(c(-0.3, 0.3))+
  # ylab("% AT")+
  ggtitle("Genes, overlapping with ELYSxNPC")

at.labs <-  seq(0.5, 0.9, by = 0.05)
at.scaled <- (at.labs * 5) - 2.4

p2 <- ggplot(elys.at %>% 
               filter(variable == "nuc"), aes(x = kb))+
  geom_line(aes(y = elys),col = "blue")+
  geom_line(aes(y = AT * 3 - 1.8), col = "red")+
  scale_y_continuous(
    name = "log2(Dam-Elys/Dam)",
    limits = c(-0.5, 0.9),
    sec.axis = sec_axis(trans = (~.),
                        breaks = at.scaled,
                        labels = at.labs,
                        name = "% AT")
  )+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "blue", size=13),
        axis.title.y.right = element_text(color = "red", size=13))+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  # ylab("% AT")+
  ggtitle("Genes, overlapping with ELYSxNUC")

at.labs <-  seq(0.4, 0.7, by = 0.05)
at.scaled <- (at.labs * 5) - 2.8

p3 <- ggplot(elys.at %>% 
               filter(variable == "other"), aes(x = kb))+
  geom_line(aes(y = elys),col = "blue")+
  geom_line(aes(y = AT * 3 - 1.8), col = "red")+
  scale_y_continuous(
    name = "log2(Dam-Elys/Dam)",
    limits = c(-0.5, 0.9),
    sec.axis = sec_axis(trans = (~.),
                        breaks = at.scaled,
                        labels = at.labs,
                        name = "% AT")
  )+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "blue", size=13),
        axis.title.y.right = element_text(color = "red", size=13))+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  # ylab("% AT")+
  ggtitle("Genes, overlapping with ELYS (no Nups)")

pdf("plots/genes.x.elys.nups.at.300bp.elys.pr.2.pdf", width = 8, height = 5)
p1
p2
p3
dev.off()

save(p1, p2, p3, file = "RData/genes.x.elys.metagene.for.elys.and.AT.plots.RData")

