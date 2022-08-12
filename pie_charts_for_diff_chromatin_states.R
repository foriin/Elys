library(reshape2)
library(stringr)
library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(openxlsx)

load("RData/elys.x.nup.genes.shev.RData", verbose = T)

# Filion 5-state color scheme
states.5 <- read.xlsx("tables/5 state chromatin model_Kc_Filion 2010.xlsx",
                      startRow = 2, colNames = F) %>% setNames(c("chr", "start", "end", "col")) %>% 
  mutate(chr = sub("chr", "", chr)) %>% filter(!is.na(end))

states.5.gr <- makeGRangesFromDataFrame(states.5, keep.extra.columns = T)

states.5.grl <- split(states.5.gr, states.5.gr$col)

findovers <- function(gr){
  sapply(states.5.grl, function(x) sum(width(GenomicRanges::intersect(gr, x, ignore.strand = T))))
}


findovers(g.npc.gr)/sum(width(g.npc.gr)) 
findovers(g.nuc.gr)/sum(width(g.nuc.gr))
findovers(g.npc.nuc.gr)/sum(width(g.npc.nuc.gr))

# Pies for Elys x npc
elys.npc.x.5states.plot <- data.frame(
  group = factor(c("Black", "Blue", "Green", "Red", "Yellow"),
                 levels = c("Black", "Blue", "Green", "Red", "Yellow")),
  value = findovers(nup.npc.x.elys.emb.1)/sum(width(nup.npc.x.elys.emb.1))
)

pie1 <- ggplot(elys.npc.x.5states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_manual(values = c("black", "blue", "green", "red", "yellow"),
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 NPC) domains\n in Filion chromatin colours")


# Pies for Elys x nuc
elys.nuc.x.5states.plot <- data.frame(
  group = factor(c("Black", "Blue", "Green", "Red", "Yellow"),
                 levels = c("Black", "Blue", "Green", "Red", "Yellow")),
  value = findovers(nup.nuc.x.elys.emb.1)/sum(width(nup.nuc.x.elys.emb.1))
)

pie2 <- ggplot(elys.nuc.x.5states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_manual(values = c("black", "blue", "green", "red", "yellow"),
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 NUC) domains\n in Filion chromatin colours")


load("RData/nup98.damid.hmm3.RData", verb = T)
nup.npc.nuc <- GenomicRanges::intersect(nup98.npc.hmm3, nup98.nuc.hmm3,
                                        ignore.strand = T)

elys.x.nups.ix <- GenomicRanges::intersect(elys.emb.hmm3.ng, nup.npc.nuc,
                                           ignore.strand = T)


# Pies for Elys x (NPC x NUC)
elys.nups.ix.x.5states.plot <- data.frame(
  group = factor(c("Black", "Blue", "Green", "Red", "Yellow"),
                 levels = c("Black", "Blue", "Green", "Red", "Yellow")),
  value = findovers(elys.x.nups.ix)/sum(width(elys.x.nups.ix))
)

pie3 <- ggplot(elys.nups.ix.x.5states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_manual(values = c("black", "blue", "green", "red", "yellow"),
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 (NPCxNUC)) domains\n in Filion chromatin colours")


pdf("plots/elys.nups.x.5.colours.pdf")
  grid.arrange(
    pie1,
    pie2,
    pie3,
    ncol = 2,
    nrow = 2
  )
dev.off()


# Kharchenko 9-state color scheme

load("RData/kharchenko.chromatin.colours.RData", verbose = T)

states.9.grl <- split(col.gr, col.gr$colour)

findovers <- function(gr){
  sapply(states.9.grl, function(x) sum(width(GenomicRanges::intersect(gr, x, ignore.strand = T))))
}


# Pies for Elys x npc
elys.npc.x.9states.plot <- data.frame(
  group = factor(1:9,
                 levels = 1:9),
  value = findovers(nup.npc.x.elys.emb.1)/sum(width(nup.npc.x.elys.emb.1))
)

pie1 <- ggplot(elys.npc.x.9states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_brewer(palette = "Set1",
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 NPC) domains\n in Kharchenko chromatin colours")

# Pies for Elys x nuc
elys.nuc.x.9states.plot <- data.frame(
  group = factor(1:9,
                 levels = 1:9),
  value = findovers(nup.nuc.x.elys.emb.1)/sum(width(nup.nuc.x.elys.emb.1))
)

pie2 <- ggplot(elys.nuc.x.9states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 1),
                                         "%")), color = "white", size=4) +
  scale_fill_brewer(palette = "Set1",
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 NUC) domains\n in Kharchenko chromatin colours")


# Pies for Elys x (NPC x NUC)
elys.nups.ix.x.9states.plot <- data.frame(
  group = factor(1:9,
                 levels = 1:9),
  value = findovers(elys.x.nups.ix)/sum(width(elys.x.nups.ix))
)

pie3 <- ggplot(elys.nups.ix.x.9states.plot %>% 
                 
                 arrange(desc(group)) %>%
                 mutate(prop = value *100) %>%
                 mutate(ypos = cumsum(prop)- 0.5*prop ),
               aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 1),
                                         "%")), color = "white", size=4) +
  scale_fill_brewer(palette = "Set1",
                    name = "Chromatin colour")+
  ggtitle("(Elys embryo x Nup98 (NPCxNUC)) domains\n in Kharchenko chromatin colours")


pdf("plots/elys.nups.x.9.colours.pdf")
grid.arrange(
  pie1,
  pie2,
  pie3,
  ncol = 2,
  nrow = 2
)
dev.off()
