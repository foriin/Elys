# 07/07/21 Here I'll use emb lam and elys domains as well as NUP* domains
# with gaps lesser than 900 nt filled

# This script is used to make averaged DamID-plots around centers
# of different genomic regions usually obtained by intersections of Elys
# domains with Nup98 domains and/or some other loci
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

load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("../Lam_sperm/RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)
load("RData/s2.elys.chip.bin.RData", verbose = T)
load("RData/h3k27ac.RData", verbose = T)

# Let's make profiles of Kc167 Lamin DamID around centers of embryonic Elys
# domains that overlap with Kc167 NPC* domains and LADs

# 1. Embryonic Elys domains that overlap w/ Kc NPC* and LADs

domains.1 <- GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                      elys.emb.hmm3.ng,
                                      ignore.strand = T)
domains.1 <- subsetByOverlaps(domains.1[width(domains.1) > 100],
                              kc.lam.hmm,
                              maxgap = 2000,
                              ignore.strand = T)
summary(width(domains.1))
# 
domains.1.c <- resize(domains.1, 1,fix = "center")

profs.1 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2,
  h3k27ac.chip.bin
)
profs.1[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.1[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.1[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)"
)
profs.1[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys"
)
profs.1[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "Emb Elys x Nup98 NPC* x Kc167 LADs\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input"
)

ps.1 <- lapply(profs.1, function(pr){
  prof  <- tss.surr(domains.1.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.3,se = T, span = 1/5, size = 1.5,
                col = "black")+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})





pdf("plots/lam.elys.profiles.around.elys.x.npc.domains.no.gap.pdf")
grid.arrange(ps.1[[1]],
             ps.1[[2]], nrow = 2)
grid.arrange(ps.1[[3]],
             ps.1[[4]], nrow = 2)
grid.arrange(ps.1[[5]], nrow = 2)
dev.off()

# 2. Embryonic Elys domains that overlap w/ Kc NPC* and embryonic LADs
# + S2 Elys Chip profiles

domains.2 <- GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                      elys.emb.hmm3.ng,
                                      ignore.strand = T) %>% 
  subset(width > 100) %>% 
  subsetByOverlaps(lam.emb.hmm2.ng,
                              maxgap = 2000,
                              ignore.strand = T)
summary(width(domains.2))

domains.2.c <- resize(domains.2, 1,fix = "center")


profs.2 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2,
  h3k27ac.chip.bin
)
profs.2[[1]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x embryo Lads\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.2[[2]]@metadata <- list(
  protein = "Lam",
  title = "Emb Elys x Nup98 NPC* x embryo Lads\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.2[[3]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x embryo Lads\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)"
)
profs.2[[4]]@metadata <- list(
  protein = "Elys",
  title = "Emb Elys x Nup98 NPC* x embryo Lads\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys"
)
profs.2[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "Emb Elys x Nup98 NPC* x embryo Lads\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input"
)

ps.2 <- lapply(profs.2, function(pr){
  prof  <- tss.surr(domains.2.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.3,se = T, span = 1/5, size = 1.5,
                col = "black")+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})


pdf("plots/lam.elys.profiles.around.elys.x.npc.domains.x.emb.lads.no.gap.pdf")
grid.arrange(ps.2[[1]],
             ps.2[[2]], nrow = 2)
grid.arrange(ps.2[[3]],
             ps.2[[4]], nrow = 2)
grid.arrange(ps.2[[5]], nrow = 2)
dev.off()


# 4. S2 Elys peaks that overlap w/ LADs in Kc 
load("RData/s2.elys.chip.RData", verbose = T)
domains.4 <- subsetByOverlaps(s2.chip.bed,
                                      kc.lam.hmm,
                                      ignore.strand = T,
                              maxgap = 2000)
summary(width(domains.4))

domains.4.c <- resize(domains.4, 1,fix = "center")


profs.4 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2,
  h3k27ac.chip.bin
)
profs.4[[1]]@metadata <- list(
  protein = "Lam",
  title = "S2 Elys x Kc167 LADs\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.4[[2]]@metadata <- list(
  protein = "Lam",
  title = "S2 Elys x Kc167 LADs\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.4[[3]]@metadata <- list(
  protein = "Elys",
  title = "S2 Elys x Kc167 LADs\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)"
)
profs.4[[4]]@metadata <- list(
  protein = "Elys",
  title = "S2 Elys x Kc167 LADs\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys"
)
profs.4[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "S2 Elys x Kc167 LADs\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input"
)

ps.4 <- lapply(profs.4, function(pr){
  prof  <- tss.surr(domains.4.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    # geom_line(aes(y = ymax), col = "red")+
    # geom_line(aes(y = ym), col = "blue")
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
      stat_smooth(geom = "line",
      method = "loess",alpha = 0.3,se = T, span = 1/5, size = 1.5,
      col = "black")+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})


pdf("plots/lam.elys.profiles.around.s2.elys.x.kc.lads.no.gap.pdf")
grid.arrange(ps.4[[1]],
             ps.4[[2]], nrow = 2)
grid.arrange(ps.4[[3]],
             ps.4[[4]], nrow = 2)
grid.arrange(ps.4[[5]], nrow = 2)
dev.off()

# Averaged profiles of Kc lam and S2 H3K27Ac around centered S2 Nup93 peaks

load("RData/s2.nups.capelson.ChIP.RData", verbose = T)

domains.5.c <- resize(s2.nup93.chip, 1, fix = "center")


profs.5 <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2,
  h3k27ac.chip.bin
)
profs.5[[1]]@metadata <- list(
  protein = "Lam",
  title = "S2 Nup93\nLam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.5[[2]]@metadata <- list(
  protein = "Lam",
  title = "S2 Nup93\nLam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs.5[[3]]@metadata <- list(
  protein = "Elys",
  title = "S2 Nup93\nElys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)"
)
profs.5[[4]]@metadata <- list(
  protein = "Elys",
  title = "S2 Nup93\nElys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys"
)
profs.5[[5]]@metadata <- list(
  protein = "H3K27Ac",
  title = "S2 Nup93\nH3K27Ac ChIP-seq in S2",
  y_lab = "ChIP-seq H3K27Ac/Input"
)

ps.5 <- lapply(profs.5, function(pr){
  prof  <- tss.surr(domains.5.c, pr, vic = 20, sig.col = "score")
  dom.prof <- prof %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    # geom_line(aes(y = ymax), col = "red")+
    # geom_line(aes(y = ym), col = "blue")
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
      method = "loess",alpha = 0.3,se = T, span = 1/5, size = 1.5,
      col = "black")+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(pr@metadata$title)
})


pdf("plots/diff.profiles.around.s2.nup93.2.pdf")
grid.arrange(ps.5[[1]],
             ps.5[[2]], nrow = 2)
grid.arrange(ps.5[[3]],
             ps.5[[4]], nrow = 2)
grid.arrange(ps.5[[5]], nrow = 2)
dev.off()

