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
load("RData/nups.x.elys.emb.RData", verbose = T)


bins <- tileGenome(seqinfo(BSgenome.Dmelanogaster.UCSC.dm3),tilewidth = 50,
                   cut.last.tile.in.chrom = T)
bins <- bins[seqnames(bins) %in% paste0("chr", euc.chroms)]
bins.seq <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3,bins)
score <- mclapply(bins.seq, function(string){
  1 - GC(s2c(as.character(string)))
}, mc.cores = 12)
seqlevels(bins) <- sub("chr", "", seqlevels(bins))

bw.pre <- bins[!is.na(bins$AT)]
export.bw(bw.pre, "bed/dm3.100bp.at.perc.bigWig")

at.npc <- read.table("bed/metaplot.elys.npc.at.count.tab", skip = 1)[, -1] %>% t()
at.nuc <- read.table("bed/metaplot.elys.nuc.at.count.tab", skip = 1)[, -1] %>% t()
at.oth <- read.table("bed/metaplot.elys.other.at.count.tab", skip = 1)[, -1] %>% t()
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

pdf("plots/metagene.AT.percent.pdf")
p.at
dev.off()


