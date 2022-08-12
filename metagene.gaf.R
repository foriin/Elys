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
library(rtracklayer)

# Metagene plots for GAF enrichment for genes, overlapping with 
# Elys x NPC, Elys x NUC or Elys without NPC or NUC domains

# load embryo GAF ChIP data
gaf.chip <- import.wig("bed/GSE180812_CHIP_GAF_interphase_700000_merge_rep1_rep2_rep3.wig") %>% 
  unlist() %>% GRanges()
dm6.dm3.chain <- import.chain("~/Database/dmel/Genome/dm6ToDm3.over.chain")
gaf.chip <- unlist(liftOver(gaf.chip, dm6.dm3.chain))
end(gaf.chip) <- start(gaf.chip) + 49
bins.300 <- tileGenome(seqlengths(BSgenome.Dmelanogaster.UCSC.dm3),
                       tilewidth = 300, cut.last.tile.in.chrom = T)
g.c.dj.1 <- disjoin(gaf.chip, with.revmap = T)
g.c.dj.2 <- g.c.dj.1[sapply(g.c.dj.1$revmap, length) < 2]
g.c.dj.2$score <- gaf.chip$score[unlist(g.c.dj.2$revmap)]
g.c.dj.2 <- keepSeqlevels(g.c.dj.2,
                          paste0("chr", names(chr.lengths)),
                          pruning.mode = 'coarse')
score <- mcolAsRleList(g.c.dj.2, "score")
gaf.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                        paste0("chr", names(chr.lengths)),
                                        pruning.mode = 'coarse'), score, "score",
                          na.rm = T)

bw.pre <- gaf.300[!is.na(gaf.300$score)]
seqlevels(bw.pre) <- sub("chr", "", seqlevels(gaf.300))
export.bw(bw.pre, "bed/emb.gaf.chip.300bp.bigWig")

# import deeptools data

gaf.npc <- read.table("metagene/metaplot.elys.npc.gaf.300.tab", skip = 1)[, -1] %>% t()
gaf.nuc <- read.table("metagene/metaplot.elys.nuc.gaf.300.tab", skip = 1)[, -1] %>% t()
gaf.oth <- read.table("metagene/metaplot.elys.other.gaf.300.tab", skip = 1)[, -1] %>% t()
gaf.pr <- as.data.frame(cbind(gaf.npc, gaf.nuc, gaf.oth))
gaf.pr <- gaf.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

gaf.m <- melt(gaf.pr, id.vars = "kb")

p.gaf <- ggplot(gaf.m, aes(x = kb, y = value))+
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
  ylab("Embryo GAF/Input")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/ (NPC x NUC)")

pdf("plots/metagene.gaf.300bp.pdf")
p.gaf
dev.off()

# load S2 GAF ChIP

s2.gaf.wig <- import.wig("bed/GAF_chip_S2_profile_GSE40646.wig")
s2.g.c.dj.1 <- disjoin(s2.gaf.wig, with.revmap = T)
s2.g.c.dj.2 <- s2.g.c.dj.1[sapply(s2.g.c.dj.1$revmap, length) < 2]
s2.g.c.dj.2$score <- s2.gaf.wig$score[unlist(s2.g.c.dj.2$revmap)]
s2.g.c.dj.2 <- keepSeqlevels(s2.g.c.dj.2,
                          paste0("chr", names(chr.lengths)),
                          pruning.mode = 'coarse')
score <- mcolAsRleList(s2.g.c.dj.2, "score")
s2.gaf.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                           paste0("chr", names(chr.lengths)),
                                           pruning.mode = 'coarse'),
                             score, "score",na.rm = T)

bw.pre <- s2.gaf.300[!is.na(s2.gaf.300$score)]
seqlevels(bw.pre) <- sub("chr", "", seqlevels(gaf.300))
export.bw(bw.pre, "bed/s2.gaf.chip.300bp.bigWig")

# import deeptools data

s2.gaf.npc <- read.table("metagene/metaplot.elys.npc.s2.gaf.300.tab", skip = 1)[, -1] %>% t()
s2.gaf.nuc <- read.table("metagene/metaplot.elys.nuc.s2.gaf.300.tab", skip = 1)[, -1] %>% t()
s2.gaf.oth <- read.table("metagene/metaplot.elys.other.s2.gaf.300.tab", skip = 1)[, -1] %>% t()
s2.gaf.pr <- as.data.frame(cbind(s2.gaf.npc, s2.gaf.nuc, s2.gaf.oth))
s2.gaf.pr <- s2.gaf.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

s2.gaf.m <- melt(s2.gaf.pr, id.vars = "kb")

p.gaf.s2 <- ggplot(s2.gaf.m, aes(x = kb, y = value))+
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
  ylab("S2 GAF ChIP")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/ (NPC x NUC)")

pdf("plots/metagene.s2.gaf.300bp.pdf")
p.gaf.s2
dev.off()

save(p.gaf, p.gaf.s2, file = "RData/elys.nups.genes.metagene.gaf.plots.RData")
