# 08/09/21 Here I'll use emb lam and elys domains,Kc NUP* and Lam domains
# with gaps lesser than 900 nt filled

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

plot.f <- function(domains, profs, ggt){
  lapply(profs, function(pr){
    prof  <- tss.surr(domains, pr, vic = 20, sig.col = "score")
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
                  method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.2,
                  col = "black")+
      ylab(pr@metadata$y_lab)+
      theme_bw()+
      scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
      ggtitle(paste0(ggt, "\n",
                     pr@metadata$title))
  })
}

load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)
load("RData/s2.elys.chip.bin.RData", verbose = T)
load("RData/cons.LADS.RData", verbose = T)

profs <- list(
  kc.bin.lam.pr,
  lam.emb.pr,
  elys.emb.pr,
  s2.elys.chip.bins.2
)
profs[[1]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs[[2]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs[[3]]@metadata <- list(
  protein = "Elys",
  title = "Elys DamID in Drosophila embryos",
  y_lab = "log2(Dam/Dam-Elys)"
)
profs[[4]]@metadata <- list(
  protein = "Elys",
  title = "Elys ChIP-seq in S2",
  y_lab = "ChIP-seq Elys"
)

# 1. Embryonic Elys domains that overlap w/ Kc NPC* and сLADs

domains.1 <- GenomicRanges::intersect(nup98.npc.hmm3.uq.ng,
                                      elys.emb.hmm3.ng,
                                      ignore.strand = T)
domains.1a <- subsetByOverlaps(domains.1[width(domains.1) > 100],
                               clads,
                               maxgap = 2000,
                               ignore.strand = T) %>% 
  resize(1,fix = "center")
domains.1b <- subsetByOverlaps(domains.1[width(domains.1) > 100],
                               clads,
                               type = 'within',
                               ignore.strand = T) %>% 
  resize(1, fix = "center")
domains.1c <- GenomicRanges::intersect(domains.1[width(domains.1) > 100],
                                       clads,
                                       ignore.strand = T) %>% 
  resize(1, fix = "center")


ps.1a <- plot.f(domains.1a, profs, "Emb Elys x Nup98 NPC* overlap cLADs +- 2kb")
ps.1b <- plot.f(domains.1b, profs, "Emb Elys x Nup98 NPC* within cLADs")
ps.1c <- plot.f(domains.1c, profs, "Emb Elys x Nup98 NPC* x cLADs")

ps.1 <- lapply(profs, function(pr){
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
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.2,
                col = "black")+
    ylab(pr@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(paste0("Emb Elys x Nup98 NPC* x cLADs\n",
                   pr@metadata$title))
})









# 2. Elys S2 ChIP-peaks that overlap w/ NPC* and cLADs

load("RData/s2.elys.chip.RData", verbose = T)
domains.2 <- GenomicRanges::intersect(s2.chip.bed,
                              nup98.npc.hmm3.uq.ng,
                              ignore.strand = T)
domains.2a <- subsetByOverlaps(domains.2[width(domains.2) > 100],
                              clads,
                              maxgap = 2000,
                              ignore.strand = T) %>% 
  resize(1,fix = "center")
domains.2b <- subsetByOverlaps(domains.2[width(domains.2) > 100],
                               clads,
                               type = 'within',
                               ignore.strand = T) %>% 
  resize(1, fix = "center")
domains.2c <- GenomicRanges::intersect(domains.2[width(domains.2) > 100],
                               clads,
                               ignore.strand = T) %>% 
  resize(1, fix = "center")


ps.2a <- plot.f(domains.2a, profs)
ps.2b <- plot.f(domains.2b, profs)
ps.2c <- plot.f(domains.2c, profs)
ps.2a <- plot.f(domains.2a, profs, "S2 Elys x Nup98 NPC* overlap cLADs +- 2kb")
ps.2b <- plot.f(domains.2b, profs, "S2 Elys x Nup98 NPC* within cLADs")
ps.2c <- plot.f(domains.2c, profs, "S2 Elys x Nup98 NPC* x cLADs")

pdf("plots/profiles.around.elys.s2.x.npc.x.clads.domains.no.gap.pdf", height = 7,
    width = 8)
grid.arrange(ps.2[[1]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black")),
             ps.2[[2]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black")), ncol = 2)
grid.arrange(ps.2[[3]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black")),
             ps.2[[4]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black")), ncol = 2)
dev.off()

pdf("plots/kc.lam.profiles.around.elys.s2.or.emb.x.npc.x.clads.domains.no.gap.pdf", height = 7,
    width = 10)
grid.arrange(ps.1a[[1]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+
               ylim(c(-0.2, 0.5))+
               annotate(
                 "text", label = paste0(length(domains.1a), " domains"),
                 x = -5, y = 0.5, size = 5
               ),
             ps.2a[[1]] + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"))+
               ylim(c(-0.2, 0.5))+
               annotate(
                 "text", label = paste0(length(domains.2a), " domains"),
                 x = -5, y = 0.5, size = 5
               ),
             ncol = 2)
grid.arrange(ps.1b[[1]] + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+
               ylim(c(-0.1, 0.65)),
             ps.2b[[1]] + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+
               ylim(c(-0.1, 0.65)),
             ncol = 2)
grid.arrange(ps.1c[[1]] + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+
               ylim(c(-0.1, 0.65)),
             ps.2c[[1]] + theme(panel.border = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.line = element_line(colour = "black"))+
               ylim(c(-0.1, 0.65)),
             ncol = 2)
dev.off()




