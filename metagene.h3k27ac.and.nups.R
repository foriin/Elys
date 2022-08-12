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

# Metagene plots for H3K27Ac enrichment for genes, overlapping with 
# Elys x NPC, Elys x NUC or Elys without NPC or NUC domains

# load S2 H3K27Ac ChIP data
load("RData/h3k27ac.RData", verbose = T)

bw.pre <- h3k27ac.chip.bin[!is.na(h3k27ac.chip.bin$score)]
export.bw(bw.pre, "bed/s2.h3k27ac.chip.300bp.bigWig")

# import deeptools data

h3k27ac.npc <- read.table("metagene/metaplot.elys.npc.h3k27ac.tab", skip = 1)[, -1] %>% t()
h3k27ac.nuc <- read.table("metagene/metaplot.elys.nuc.h3k27ac.tab", skip = 1)[, -1] %>% t()
h3k27ac.oth <- read.table("metagene/metaplot.elys.other.h3k27ac.tab", skip = 1)[, -1] %>% t()
h3k27ac.pr <- as.data.frame(cbind(h3k27ac.npc, h3k27ac.nuc, h3k27ac.oth))
h3k27ac.pr <- h3k27ac.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

h3k27ac.m <- melt(h3k27ac.pr, id.vars = "kb")

p.h3k27ac <- ggplot(h3k27ac.m, aes(x = kb, y = value))+
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
  ylab("H3K27Ac/Input")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/ (NPC x NUC)")

pdf("plots/metagene.h3k27ac.percent.300bp.pdf")
p.h3k27ac
dev.off()

save(p.h3k27ac, file = "RData/metagene.elys.nup.x.genes.h3k27ac.plot.RData")

# Metagene plots for nup profiles?
load("RData/nup98.damid.log2.RData", verb = T)
nup98.npc.wig
# Nup98 NPC
seqlevels(nup98.npc.wig) <- sub(" ", "", seqlevels(nup98.npc.wig))
nup.nuc.score <- mcolAsRleList(nup98.nucl.wig, "score")
bins.300 <- tileGenome(seqlengths(BSgenome.Dmelanogaster.UCSC.dm3),
                       tilewidth = 300, cut.last.tile.in.chrom = T)
nup.npc.dj.1 <- disjoin(nup98.npc.wig, with.revmap = T)
nup.npc.dj.2 <- nup.npc.dj.1[sapply(nup.npc.dj.1$revmap, length) < 2]
nup.npc.dj.2$score <- nup98.npc.wig$score[unlist(nup.npc.dj.2$revmap)]
nup.npc.dj.2 <- keepSeqlevels(nup.npc.dj.2,
                           paste0("chr", names(chr.lengths)),
                           pruning.mode = 'coarse')
score.nup.npc <- mcolAsRleList(nup.npc.dj.2, "score")
nup.npc.damid.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                             seqlevels(nup.npc.dj.2),
                                               pruning.mode = 'coarse'), score.nup.npc, "score",
                                 na.rm = T)
bw.pre.npc <- nup.npc.damid.300[!(is.na(nup.npc.damid.300$score))]
export.bw(bw.pre.npc, "bed/kc.nup98.npc.damid.300bp.bigWig")


# import deeptools data

nup98.npc.npc <- read.table("metagene/nup98/metaplot.elys.npc.nup98.npc.300.tab", skip = 1)[, -1] %>% t()
nup98.npc.nuc <- read.table("metagene/nup98/metaplot.elys.nuc.nup98.npc.300.tab", skip = 1)[, -1] %>% t()
nup98.npc.oth <- read.table("metagene/nup98/metaplot.elys.other.nup98.npc.300.tab", skip = 1)[, -1] %>% t()
nup98.npc.pr <- as.data.frame(cbind(nup98.npc.npc, nup98.npc.nuc, nup98.npc.oth))
nup98.npc.pr <- nup98.npc.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

nup98.npc.m <- melt(nup98.npc.pr, id.vars = "kb")

p.nup98.npc <- ggplot(nup98.npc.m, aes(x = kb, y = value))+
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
  ylab("log2(Nup98_NPC/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/ (NPC x NUC)")

pdf("plots/metaplot.elys.x.nups.nup98.npc.300bp.pdf")
p.nup98.npc
dev.off()

save(p.nup98.npc, file = "RData/elys.nup.x.genes.nup98.npc.plot.RData")


# Nup98 Nucleoplasmic
seqlevels(nup98.nucl.wig) <- sub(" ", "", seqlevels(nup98.nucl.wig))
bins.300 <- tileGenome(seqlengths(BSgenome.Dmelanogaster.UCSC.dm3),
                       tilewidth = 300, cut.last.tile.in.chrom = T)
nup.nucl.dj.1 <- disjoin(nup98.nucl.wig, with.revmap = T)
nup.nucl.dj.2 <- nup.nucl.dj.1[sapply(nup.nucl.dj.1$revmap, length) < 2]
nup.nucl.dj.2$score <- nup98.nucl.wig$score[unlist(nup.nucl.dj.2$revmap)]
nup.nucl.dj.2 <- keepSeqlevels(nup.nucl.dj.2,
                              paste0("chr", names(chr.lengths)),
                              pruning.mode = 'coarse')
score.nup.nucl <- mcolAsRleList(nup.nucl.dj.2, "score")
nup.nucl.damid.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                                  seqlevels(nup.nucl.dj.2),
                                                  pruning.mode = 'coarse'), score.nup.nucl, "score",
                                    na.rm = T)
bw.pre.nucl <- nup.nucl.damid.300[!(is.na(nup.nucl.damid.300$score))]
export.bw(bw.pre.nucl, "bed/kc.nup98.nucl.damid.300bp.bigWig")

nup98.nuc.npc <- read.table("metagene/nup98/metaplot.elys.npc.nup98.nuc.300.tab", skip = 1)[, -1] %>% t()
nup98.nuc.nuc <- read.table("metagene/nup98/metaplot.elys.nuc.nup98.nuc.300.tab", skip = 1)[, -1] %>% t()
nup98.nuc.oth <- read.table("metagene/nup98/metaplot.elys.other.nup98.nuc.300.tab", skip = 1)[, -1] %>% t()
nup98.nuc.pr <- as.data.frame(cbind(nup98.nuc.npc, nup98.nuc.nuc, nup98.nuc.oth))
nup98.nuc.pr <- nup98.nuc.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

nup98.nuc.m <- melt(nup98.nuc.pr, id.vars = "kb")

p.nup98.nuc <- ggplot(nup98.nuc.m, aes(x = kb, y = value))+
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
  ylab("log2(Nup98_nuc/Dam)")+
  ggtitle("Genes, overlapping with ELYSxnuc, ELYSxNUC or ELYS w/ (nuc x NUC)")

library(gridExtra)
pdf("plots/metaplot.elys.x.nups.nup98.damid.300bp.pdf", height = 8)
grid.arrange(
  p.nup98.npc,
  p.nup98.nuc,
  nrow = 2
)
dev.off()

save(p.nup98.nuc, file = "RData/elys.nup.x.genes.nup98.nuc.plot.RData")


