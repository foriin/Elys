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

# Metagene plots for Pol II enrichment for genes, overlapping with 
# Elys x NPC, Elys x NUC or Elys without NPC or NUC domains

# load S2 PolII ChIP data

polII.chip <- import.bw("bed/GSM4745381_ChIP-Pol-S2P_DMSO_mock-RNAi_S2.bigwig")
polII.input <- import.bw("bed/GSM4745389_Input_DMSO_mock-RNAi_S2.bigwig")
dm6.dm3.chain <- import.chain("~/Database/dmel/Genome/dm6ToDm3.over.chain")
polII.chip.dm3 <- unlist(liftOver(polII.chip, dm6.dm3.chain))
polII.input.dm3 <- unlist(liftOver(polII.input, dm6.dm3.chain))
bins.300 <- tileGenome(seqlengths(BSgenome.Dmelanogaster.UCSC.dm3),
                       tilewidth = 300, cut.last.tile.in.chrom = T)
pol2.dj.1 <- disjoin(polII.chip.dm3, with.revmap = T)
pol2.dj.2 <- pol2.dj.1[sapply(pol2.dj.1$revmap, length) < 2]
pol2.dj.2$score <- polII.chip.dm3$score[unlist(pol2.dj.2$revmap)]
pol2.dj.2 <- keepSeqlevels(pol2.dj.2,
                          paste0("chr", names(chr.lengths)),
                          pruning.mode = 'coarse')
score.chip <- mcolAsRleList(pol2.dj.2, "score")
polII.chip.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                        paste0("chr", names(chr.lengths)),
                                        pruning.mode = 'coarse'), score.chip, "score",
                          na.rm = T)
input.dj.1 <- disjoin(polII.input.dm3, with.revmap = T)
input.dj.2 <- input.dj.1[sapply(input.dj.1$revmap, length) < 2]
input.dj.2$score <- polII.input.dm3$score[unlist(input.dj.2$revmap)]
input.dj.2 <- keepSeqlevels(input.dj.2,
                           paste0("chr", names(chr.lengths)),
                           pruning.mode = 'coarse')
score.input <- mcolAsRleList(input.dj.2, "score")
polII.input.300 <-  binnedAverage(keepSeqlevels(bins.300,
                                               paste0("chr", names(chr.lengths)),
                                               pruning.mode = 'coarse'), score.input, "score",
                                 na.rm = T)


pol2.chip.300 <- polII.chip.300
pol2.chip.300$score <- polII.chip.300$score/(polII.input.300$score + 1)

bw.pre <- pol2.chip.300[!is.na(pol2.chip.300$score)]
seqlevels(bw.pre) <- sub("chr", "", seqlevels(bw.pre))
export.bw(bw.pre, "bed/s2.pol2.chip.300bp.bigWig")

# import deeptools data

pol2.npc <- read.table("metagene/metaplot.elys.npc.pol2.300.tab", skip = 1)[, -1] %>% t()
pol2.nuc <- read.table("metagene/metaplot.elys.nuc.pol2.300.tab", skip = 1)[, -1] %>% t()
pol2.oth <- read.table("metagene/metaplot.elys.other.pol2.300.tab", skip = 1)[, -1] %>% t()
pol2.pr <- as.data.frame(cbind(pol2.npc, pol2.nuc, pol2.oth))
pol2.pr <- pol2.pr %>% setNames(c("npc", "nuc", "other")) %>% 
  mutate(kb = 1:600)

pol2.m <- melt(pol2.pr, id.vars = "kb")

p.pol2 <- ggplot(pol2.m, aes(x = kb, y = value))+
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
  ylab("PolII/Input")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC or ELYS w/ (NPC x NUC)")

pdf("plots/metagene.pol2.chip.s2.300bp.pdf")
p.pol2
dev.off()

save(p.pol2, file = "RData/metagene.elys.nup.x.genes.polII.plot.RData")
