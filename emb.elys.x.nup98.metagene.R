library(dm3)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(tibble)
library(eulerr)
library(reshape2)
library(openxlsx)

load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/emb.elys.hmm3.x.nup98.npc.nuc.uq.RData", verb = T)

# Metagene plots

# Prepare and export data for Deeptools software
export.bed(subsetByOverlaps(dm3.genes.gr, elys.x.nup.npc), "bed/dm3.genes.elys.npc.bed")
export.bed(subsetByOverlaps(dm3.genes.gr, elys.x.nup.nuc), "bed/dm3.genes.elys.nuc.bed")
dm3.rand.gr <- bedTools.shuffle.gr(dm3.genes.gr[seqnames(dm3.genes.gr) %in% euc.chroms])
export.bed(subsetByOverlaps(dm3.rand.gr, elys.emb.hmm3.ng, ignore.strand = T),
           "bed/rand.genes.x.elys.bed")
bw.pre <- elys.emb.pr[!is.na(elys.emb.pr$score)]
seqlengths(bw.pre) <- chr.lengths[match(names(seqlengths(bw.pre)), names(chr.lengths))]
export.bw(bw.pre, "bed/elys.embryos.bigWig")
bw.pre.lam <- lam.emb.pr[!is.na(lam.emb.pr$score)]
seqlengths(bw.pre.lam) <- chr.lengths[match(names(seqlengths(bw.pre.lam)), names(chr.lengths))]
export.bw(bw.pre.lam, "bed/lam.embryos.bigWig")
# Load data from deeptools computematrix and plotprofile pipeline
lam.nuc <- read.table("bed/metaplot.lam.elys.nuc.tab", skip = 1)[, -1] %>% t()
lam.npc <- read.table("bed/metaplot.lam.elys.npc.tab", skip = 1)[, -1] %>% t()
elys.nuc <- read.table("bed/metaplot.elys.nuc.tab", skip = 1)[, -1] %>% t()
elys.npc <- read.table("bed/metaplot.elys.npc.tab", skip = 1)[, -1] %>% t()


npc.pr <- as.data.frame(cbind(elys.npc, lam.npc))
npc.pr <- npc.pr %>% setNames(c("elys", "lamin")) %>% 
  mutate(kb = 1:600)
nuc.pr <- as.data.frame(cbind(elys.nuc, lam.nuc))
nuc.pr <- nuc.pr %>% setNames(c("elys", "lamin")) %>% 
  mutate(kb = 1:600)
elys.pr <- as.data.frame(cbind(elys.nuc, elys.npc))
elys.pr <- elys.pr %>% setNames(c("nuc", "npc")) %>% 
  mutate(kb = 1:600)

npc.m <- melt(npc.pr, id.vars = "kb")
nuc.m <- melt(nuc.pr, id.vars = "kb")
elys.m <- melt(elys.pr, id.vars = "kb")

pdf("plots/genes.elys.x.npc.pdf")
ggplot(npc.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("purple",
                                "red"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("log2(Dam-X/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNPC")
dev.off()

pdf("plots/genes.elys.x.nuc.pdf")
ggplot(nuc.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("purple",
                                "red"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("log2(Dam-X/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNUC")
dev.off()

pdf("plots/genes.elys.x.npc.or.nuc.pdf")
ggplot(elys.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("blue",
                                "green"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("log2(Dam-Elys/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNPC or ELYSxNUC")
dev.off()


# export.bed(nup.npc.x.elys.emb.1, "bed/elys.x.nup.npc.bed")
# export.bed(nup.nuc.x.elys.emb.1, "bed/elys.x.nup.nuc.bed")

# Here I used gene lists provided by YY

genes.npc <- read.xlsx("tables/Genes overlapped with only Elys_NPC or Elys_nucl.xlsx", sheet = 1) %>% 
  dplyr::slice(1:(n() - 1))
genes.nuc <- read.xlsx("tables/Genes overlapped with only Elys_NPC or Elys_nucl.xlsx", sheet = 2) %>% 
  dplyr::slice(1:(n() - 1))
genes.elys <- read.xlsx("tables/Genes overlapped with only Elys_NPC or Elys_nucl and with their combinations.xlsx",
                        sheet = 3) %>% 
  dplyr::slice(1:(n() - 1))


g.npc.gr <- makeGRangesFromDataFrame(genes.npc)
g.nuc.gr <- makeGRangesFromDataFrame(genes.nuc)
g.elys.gr <- makeGRangesFromDataFrame(genes.elys)
# g.npc.nuc.gr <- g.elys.gr
# save(g.npc.gr, g.nuc.gr, g.npc.nuc.gr,
#      file = "RData/elys.x.nup.genes.shev.RData")

export.bed(g.npc.gr, "bed/genes.elys.npc.bed")
export.bed(g.nuc.gr, "bed/genes.elys.nuc.bed")
export.bed(g.elys.gr, "bed/genes.elys.bed")

# Load deeptools data
lam.nuc <- read.table("metagene/metaplot.lam.elys.nuc.tab", skip = 1)[, -1] %>% t()
lam.npc <- read.table("metagene/metaplot.lam.elys.npc.tab", skip = 1)[, -1] %>% t()
elys.nuc <- read.table("metagene/metaplot.elys.nuc.tab", skip = 1)[, -1] %>% t()
elys.npc <- read.table("metagene/metaplot.elys.npc.tab", skip = 1)[, -1] %>% t()
elys.elys <- read.table("metagene/metaplot.elys.tab", skip = 1)[, -1] %>% t()
elys.rand <- read.table("metagene/metaplot.elys.random.emb.elys.tab", skip = 1)[, -1] %>% t()
lam.rand <- read.table("metagene/metaplot.elys.random.emb.lam.tab", skip = 1)[, -1] %>% t()
rand.pr <- as.data.frame(cbind(elys.rand, lam.rand))
rand.pr <- rand.pr %>% setNames(c("Elys", "Lamin")) %>% 
  mutate(kb = 1:600)
rand.m <- melt(rand.pr, id.vars = "kb")

npc.pr <- as.data.frame(cbind(elys.npc, lam.npc))
npc.pr <- npc.pr %>% setNames(c("elys", "lamin")) %>% 
  mutate(kb = 1:600)
nuc.pr <- as.data.frame(cbind(elys.nuc, lam.nuc))
nuc.pr <- nuc.pr %>% setNames(c("elys", "lamin")) %>% 
  mutate(kb = 1:600)
elys.pr <- as.data.frame(cbind(elys.nuc, elys.npc, elys.elys, elys.rand))
elys.pr <- elys.pr %>% setNames(c("nuc", "npc", "other", "random")) %>% 
  mutate(kb = 1:600)
lam.pr <- as.data.frame(cbind(lam.nuc, lam.npc))
lam.pr <- lam.pr %>% setNames(c("nuc", "npc")) %>% 
  mutate(kb = 1:600)
npc.m <- melt(npc.pr, id.vars = "kb")
nuc.m <- melt(nuc.pr, id.vars = "kb")
elys.m <- melt(elys.pr, id.vars = "kb")
lam.m <- melt(lam.pr, id.vars = "kb")

p.lam <- ggplot(lam.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("purple",
                                "red"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("log2(Dam-Lam/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNUC or ELYSxNPC")
pdf("plots/genes.elys.x.npc.or.nuc.Lam.pr.pdf")
p.lam
dev.off()

p.elys <- ggplot(elys.m, aes(x = kb, y = value))+
  geom_line(aes(col = variable))+
  scale_color_manual(values = c("blue",
                                "green",
                                "purple",
                                "darkorange"))+
  theme_bw()+
  theme(axis.title.x=element_blank(),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_x_continuous(breaks=c(1,50,550,600),
                     labels=c("-500", "start", "end", "+500"))+
  ylab("log2(Dam-Elys/Dam)")+
  ggtitle("Genes, overlapping with ELYSxNPC, ELYSxNUC, or smth other")
pdf("plots/genes.elys.x.npc.or.nuc.Elys.pr.pdf")
p.elys
dev.off()

save(g.npc.gr, g.nuc.gr, g.elys.gr, file = "RData/metagene.log2.genes.gr.RData")
save(p.elys, p.lam, file = "RData/metagene.log2.plots.RData")
save(elys.pr, lam.pr, file = "RData/elys.nup.npc.metagene.RData")

# Pie charts for 5', 3', and gene bodies

# Embryonic Elys x Nup98 NPC

g.npc.gr.5 <- resize(g.npc.gr, fix = "start", width = 1)
g.npc.gr.5$id <- 1:length(g.npc.gr)
g.npc.gr.3 <- resize(g.npc.gr, fix = "end", width = 1)
g.npc.gr.3$id <- 1:length(g.npc.gr)


a1 <- subsetByOverlaps(g.npc.gr.5, nup.npc.x.elys.emb.1, ignore.strand = T)
a2 <- subsetByOverlaps(g.npc.gr.3, nup.npc.x.elys.emb.1, ignore.strand = T)
a3 <- base::setdiff(g.npc.gr.5$id, c(a1$id,a2$id))
throw <- unique(c(a1$id, a2$id)[duplicated(c(a1$id, a2$id))])
a1 <- a1[!(a1$id %in% throw)]
a2 <- a2[!(a2$id %in% throw)]


elys.npc.x.genes.plot <- data.frame(
  group = factor(c("Only at 5'", "Only at 3'", "At both ends", "Other"),
                 levels = c("Only at 5'", "Only at 3'", "At both ends", "Other")),
  value = c(length(a1)/length(g.npc.gr),
            length(a2)/length(g.npc.gr),
            length(c(throw))/length(g.npc.gr),
            length(a3)/length(g.npc.gr))
)

pie1 <- ggplot(elys.npc.x.genes.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_brewer(palette="Set1", name = "Elys x NPC relative to genes it overlaps")+
  ggtitle("Genes overlapping (Elys embryo x Nup98 NPC)")

# Embryonic Elys x Nup98 Nucleoplasmic

g.nuc.gr.5 <- resize(g.nuc.gr, fix = "start", width = 1)
g.nuc.gr.5$id <- 1:length(g.nuc.gr)
g.nuc.gr.3 <- resize(g.nuc.gr, fix = "end", width = 1)
g.nuc.gr.3$id <- 1:length(g.nuc.gr)


a1 <- subsetByOverlaps(g.nuc.gr.5, nup.nuc.x.elys.emb.1, ignore.strand = T)
a2 <- subsetByOverlaps(g.nuc.gr.3, nup.nuc.x.elys.emb.1, ignore.strand = T)
a3 <- base::setdiff(g.nuc.gr.5$id, c(a1$id,a2$id))
throw <- unique(c(a1$id, a2$id)[duplicated(c(a1$id, a2$id))])
a1 <- a1[!(a1$id %in% throw)]
a2 <- a2[!(a2$id %in% throw)]


elys.nuc.x.genes.plot <- data.frame(
  group = factor(c("Only at 5'", "Only at 3'", "At both ends", "Other"),
                 levels = c("Only at 5'", "Only at 3'", "At both ends", "Other")),
  value = c(length(a1)/length(g.nuc.gr),
            length(a2)/length(g.nuc.gr),
            length(c(throw))/length(g.nuc.gr),
            length(a3)/length(g.nuc.gr))
)

pie2 <- ggplot(elys.nuc.x.genes.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=5) +
  scale_fill_brewer(palette="Set1", name = "Elys x NUC relative to genes it overlaps")+
  ggtitle("Genes overlapping (Elys embryo x Nup98 NUC)")

pdf("plots/elys.x.nups.relat.to.genes.pdf", width = 8, height = 5)
pie1
pie2
dev.off()

save(pie1, pie2, file = "RData/elys.x.nups.genes.pie.charts.RData")
