library(dm3)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(gridExtra)

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

plot.f <- function(domains, prof, ggt){
  prof.surr  <- tss.surr(domains, prof, vic = 20, sig.col = "score")
  dom.prof <- prof.surr %>% rowMeans.2()  %>% 
    as.data.frame()  %>% 
    mutate(kb = seq(-6, by = 0.3, length.out = 41))
  dom.prof <- cbind(dom.prof, 
                    apply(prof.surr, 1, mean.conf) %>% t %>% 
                      as.data.frame() )
  
  
  ggplot(dom.prof %>% 
           setNames(c("score", "kb", "ymin", "ymax")),
         aes(x = kb, y = score))+
    # geom_line()+
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4)+
    stat_smooth(geom = "line",
                method = "loess",alpha = 0.9,se = T, span = 1/5, size = 1.2,
                col = "black")+
    ylab(prof@metadata$y_lab)+
    theme_bw()+
    scale_x_continuous(breaks = seq(-6, 6, 1), expand = c(0.01,0.01))+
    ggtitle(paste0(ggt, "\n",
                   prof@metadata$title))

}

load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
load("RData/br.fb.glia.nrn.lads.RData", verbose = T)
load("RData/br.fb.neur.glia.lam.profiles.RData", verbose = T)
load("RData/kc.lam.profile.ma.RData", verbose = T)
load("RData/kharchenko.chromatin.colours.RData", verbose = T)

profs <- mget(ls(pattern = "pr.gr"))
names(profs)

profs[[1]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in whole brains",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs[[2]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in fat bodies",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs[[3]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in glia",
  y_lab = "log2(Dam/Dam-Lam)"
)
profs[[4]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in neurons",
  y_lab = "log2(Dam/Dam-Lam)"
)

npc.elys.lad.br <- subsetByOverlaps(nup.npc.x.elys.emb.1[width(nup.npc.x.elys.emb.1) > 100],
                                    br.lam.hmm, ignore.strand = T,
                                    maxgap = 1999) %>% 
  resize(width = 1, fix = 'center')

ps.1 <- plot.f(npc.elys.lad.br, profs[[1]], "(Emb Elys x Nup98 NPC) x (Brain LADs +- 2kb)")

npc.elys.lad.fb <- subsetByOverlaps(nup.npc.x.elys.emb.1[width(nup.npc.x.elys.emb.1) > 100],
                                    fb.lam.hmm, ignore.strand = T,
                                    maxgap = 1999) %>% 
  resize(width = 1, fix = 'center')

ps.2 <- plot.f(npc.elys.lad.fb, profs[[2]], "(Emb Elys x Nup98 NPC) x (Fat Bodies LADs +- 2kb)")

npc.elys.lad.glia <- subsetByOverlaps(nup.npc.x.elys.emb.1[width(nup.npc.x.elys.emb.1) > 100],
                                    glia.lam.hmm, ignore.strand = T,
                                    maxgap = 1999) %>% 
  resize(width = 1, fix = 'center')

ps.3 <- plot.f(npc.elys.lad.glia, profs[[3]], "(Emb Elys x Nup98 NPC) x (Glial LADs +- 2kb)")

npc.elys.lad.nrn <- subsetByOverlaps(nup.npc.x.elys.emb.1[width(nup.npc.x.elys.emb.1) > 100],
                                      nrn.lam.hmm, ignore.strand = T,
                                      maxgap = 1999) %>% 
  resize(width = 1, fix = 'center')

ps.4 <- plot.f(npc.elys.lad.nrn, profs[[4]], "(Emb Elys x Nup98 NPC) x (Neuronal LADs +- 2kb)")

pdf("plots/lam.profiles.around.emb.elys.x.nup.npc.x.lads.pdf", height = 7,
    width = 8)
grid.arrange(ps.1 + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)),
             ps.2 + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)), ncol = 2)
grid.arrange(ps.3 + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)),
             ps.4 + theme(panel.border = element_blank(),
                               panel.grid.major = element_blank(),
                               panel.grid.minor = element_blank(),
                               axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)), ncol = 2)
dev.off()


# Lam profiles in Kc LADs random spots

kc.lam.sn <- tile(kc.lam.hmm, width = 1) %>% unlist()

kc.lad.rand <- kc.lam.sn[sample(1:length(kc.lam.sn), length(npc.elys.lad.kc))]

profs[[5]] <- kc.bins.lam.pr.f
profs[[5]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in Kc167",
  y_lab = "log2(Dam/Dam-Lam)"
)

ps.5 <- plot.f(kc.lad.rand, profs[[5]], "Random positions in Kc LADs")

pdf("plots/kc.lam.pr.around.random.lad.spots.pdf", height = 7, width = 4)
ps.5 + theme(panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            plot.title = element_text(size = 10))
dev.off()

# Emb Elys x NUP NPC x s2 active chromatin

npc.elys.act.c<- subsetByOverlaps(nup.npc.x.elys.emb.1[width(nup.npc.x.elys.emb.1) > 100],
                                    col.act.st, ignore.strand = T,
                                    type = 'within') %>% resize(width = 1,
                                                                fix = 'center')

ps.6 <- plot.f(npc.elys.act.c, profs[[5]], "(Emb Elys x Nup98 NPC) x S2 active")

pdf("plots/kc.lam.pr.around.elys.npc.s2.active.states.pdf", height = 7, width = 4)
ps.6 + theme(panel.border = element_blank(),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.line = element_line(colour = "black"),
             plot.title = element_text(size = 10))
dev.off()



ps.7 <- plot.f(npc.elys.act.c, profs[[1]], "(Emb Elys x Nup98 NPC) x S2 active")
ps.8 <- plot.f(npc.elys.act.c, profs[[2]], "(Emb Elys x Nup98 NPC) x S2 active")
ps.9 <- plot.f(npc.elys.act.c, profs[[3]], "(Emb Elys x Nup98 NPC) x S2 active")
ps.10 <- plot.f(npc.elys.act.c, profs[[4]], "(Emb Elys x Nup98 NPC) x S2 active")

pdf("plots/lam.profiles.around.emb.elys.x.nup.npc.x.s2.act.pdf", height = 7,
    width = 8)
grid.arrange(ps.7 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)),
             ps.8 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)), ncol = 2)
grid.arrange(ps.9 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)),
             ps.10 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)), ncol = 2)
dev.off()


# Genes with high expression in neurons

load("RData/glia.neurons.TPMs.RData", verbose = T)
neur.genes <- neur.tpm.gr[neur.tpm.gr$TPM > 10]
npc.elys.n.g <- subsetByOverlaps(nup.npc.x.elys.emb.1,
                                 neur.genes, ignore.strand = T) %>% 
  resize(width = 1, fix = 'center')
ps.11 <- plot.f(npc.elys.n.g, profs[[4]], "(Emb Elys x Nup98 NPC) x Neuronal expressed genes")

glia.genes <- glia.tpm.gr[glia.tpm.gr$TPM > 10]
npc.elys.g.g <- subsetByOverlaps(nup.npc.x.elys.emb.1,
                                 glia.genes, ignore.strand = T) %>% 
  resize(width = 1, fix = 'center')

ps.12 <- plot.f(npc.elys.g.g, profs[[3]], "(Emb Elys x Nup98 NPC) x Glial expressed genes")

pdf("plots/lam.profiles.around.emb.elys.x.nup.npc.x.neur.glia.exp.genes.pdf", height = 7,
    width = 8)
grid.arrange(ps.11 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)),
             ps.12 + theme(panel.border = element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          axis.line = element_line(colour = "black"),
                          plot.title = element_text(size = 10)), ncol = 2)
dev.off()


save(profs, file = "RData/lam.profiles.in.br.nrn.gl.fb.kc.RData")

