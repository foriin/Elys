library(GenomicRanges)
library(dplyr)
library(dm3)
library(rtracklayer)
library(ggplot2)
library(gridExtra)
library(openxlsx)

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

jacc.ov <- function(gr1,gr2){
  ov <- findOverlaps(gr1,gr2)
  jac = (unique(ov@from) %>% length + unique(ov@to) %>% length)/(length(gr1) + length(gr2))
  jac
}

perm.test.ov <- function(gr.1, gr.2, shuffle = "ALL",
                         ncores = 14, N = 10000,
                         genome = "/home/artem/Database/dmel/Genome/dm3.genome"){
  cval.1 <- length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T))
  cval.2 <- length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T))
  ptest <- mclapply(1:N, function(somebodyoncetoldmetheworldsgonnarollme){
    if (shuffle == "ALL") {
      gr.1 <- bedTools.shuffle.gr(gr.1, genome = genome)
      gr.2 <- bedTools.shuffle.gr(gr.2, genome = genome)
    }
    else if (shuffle %in% c("1", "2")) {
      assign(paste0("gr.", shuffle),
             bedTools.shuffle.gr(get(paste0("gr.", shuffle))),
             genome = genome)
    }
    c(length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
      length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)) > cval.1,
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)) > cval.2)
  }, mc.cores = ncores)
  ptest <- do.call(rbind, ptest)
  pv.1 <- sum(ptest[,3])/N
  pv.2 <- sum(ptest[,4])/N
  avg.rand.len.1 <- mean(ptest[,1])
  avg.rand.len.2 <- mean(ptest[,2])
  c("AxB" = cval.1,
    "AxB random" = avg.rand.len.1,
    "%A" = cval.1/length(gr.1)*100,
    "BxA" = cval.2,
    "BxA random" = avg.rand.len.2,
    "%B" = cval.2/length(gr.2)*100,
    "Pval A" = pv.1,
    "Pval B" = pv.2)
}
perm.test.jac.ov <- function(gr.1, gr.2, shuffle = "ALL",
                         ncores = 14, N = 10000,
                         genome = "/home/artem/Database/dmel/Genome/dm3.genome"){
  cval <- jacc.ov(gr.1, gr.2)
  ptest <- mclapply(1:N, function(somebodyoncetoldmetheworldsgonnarollme){
    if (shuffle == "ALL") {
      gr.1 <- bedTools.shuffle.gr(gr.1, genome = genome)
      gr.2 <- bedTools.shuffle.gr(gr.2, genome = genome)
    }
    else if (shuffle %in% c("1", "2")) {
      assign(paste0("gr.", shuffle),
             bedTools.shuffle.gr(get(paste0("gr.", shuffle))),
             genome = genome)
    }
    c(length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
      length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
      jacc.ov(gr.1, gr.2) > cval)
  }, mc.cores = ncores)
  ptest <- do.call(rbind, ptest)
  pv.1 <- sum(ptest[,3])/N
  avg.rand.len.1 <- mean(ptest[,1])
  avg.rand.len.2 <- mean(ptest[,2])
  c("AxB" = length(subsetByOverlaps(gr.1, gr.2, ignore.strand = T)),
    "AxB random" = avg.rand.len.1,
    "%A" = length(subsetByOverlaps(gr.1, gr.2,
                                   ignore.strand = T))/length(gr.1)*100,
    "BxA" = length(subsetByOverlaps(gr.2, gr.1, ignore.strand = T)),
    "BxA random" = avg.rand.len.2,
    "%B" =length(subsetByOverlaps(gr.2, gr.1,
                                  ignore.strand = T))/length(gr.2)*100,
    "Pval" = pv.1)
}

load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("../Lam_sperm/RData/kc.and.nrn.lam.hmm.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)
load("RData/neur.enhancers.RData", verbose = T)
load("RData/starr-seq.ecd.enhancers.RData", verbose = T)
load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/brain.enhancers.RData", verb = T)

# Use only euchromatic genes
dm3.tss.gr <- dm3.tss.gr %>% delete.het.gr()
dm3.genes.gr <- dm3.genes.gr %>% delete.het.gr()
# Script used to find intersections of starr ecdysone-indused enhancers in S2 or 
# brain enhancers with Elys NPC* or NUC* domains

# 1. Ecdysone starr enhancers

e.npc.starr.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, starr.ecd)
# names(e.npc.starr.pv) <- c("AxB", "AxB random", "%A", "BxA", "BxA random",
                           # "%B", "Pval A", "Pval B")
# starr.npc.pv <- perm.test.ov(starr.ecd, nup98.npc.hmm3.uq.ng)

e.nuc.starr.pv <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, starr.ecd)
# names(e.nuc.starr.pv) <- c("AxB", "AxB random", "%A", "BxA", "BxA random",
                           # "%B", "Pval A", "Pval B")
# starr.nuc.pv <- perm.test.ov(starr.ecd, nup98.nuc.hmm3.uq.ng)

# 2. Brain enhancers

sum(width(col.inact.st))
s2.inact.genome <- tapply(col.inact.st, seqnames(col.inact.st), function(x) sum(width(x)))
cat(paste(names(s2.inact.genome), s2.inact.genome, sep = "\t"), sep = "\n",
    file = "bed/s2.inact.state.chrome.size")
s2.act.genome <- tapply(col.act.st, seqnames(col.act.st), function(x) sum(width(x)))
cat(paste(names(s2.act.genome), s2.act.genome, sep = "\t"), sep = "\n",
    file = "bed/s2.act.state.chrome.size")

e.npc.b.e.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, brain.enh)
e.npc.b.e.ic.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, brain.enh,
                                 genome = "/home/artem/R/projects/Elys/bed/s2.inact.state.chrome.size")
e.nuc.b.e.pv <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, brain.enh)
# b.e.npc.pv <- perm.test.ov(brain.enh, nup98.npc.hmm3.uq.ng)
# b.e.nuc.pv <- perm.test.ov(brain.enh, nup98.nuc.hmm3.uq.ng)

# 3. TSSs, gene bodies and intergenic regions
e.npc.x.tss <- GenomicRanges::intersect(nup.npc.x.elys.emb.1, dm3.tss.gr, ignore.strand = T)
e.nuc.x.tss <- GenomicRanges::intersect(nup.nuc.x.elys.emb.1, dm3.tss.gr,ignore.strand = T)

e.npc.tss.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, dm3.tss.gr)
e.nuc.tss.pv <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, dm3.tss.gr)

# e.nuc.tss.pv.2 <- perm.test.ov(dm3.tss.gr, nup.nuc.x.elys.emb.1)
# e.npc.tss.pv.2 <- perm.test.ov(dm3.tss.gr, nup.npc.x.elys.emb.1)

e.npc.dom.1.b <- subsetByOverlaps(nup.npc.x.elys.emb.1,
                                  e.npc.x.tss, ignore.strand = T, invert = T)
e.nuc.dom.1.b <- subsetByOverlaps(nup.nuc.x.elys.emb.1,
                                  e.nuc.x.tss, ignore.strand = T, invert = T)
genes.subs.npc <- subsetByOverlaps(dm3.genes.gr, e.npc.x.tss,invert = T,
                                   ignore.strand = T)
genes.subs.nuc <- subsetByOverlaps(dm3.genes.gr,
                                   e.nuc.x.tss, invert = T, ignore.strand = T)
e.npc.gb.pv <- perm.test.jac.ov(e.npc.dom.1.b, genes.subs.npc)
e.nuc.gb.pv <- perm.test.jac.ov(e.nuc.dom.1.b, genes.subs.nuc)

intergenes <- gaps(dm3.genes.gr)
e.npc.x.gb <- GenomicRanges::intersect(e.npc.dom.1.b, genes.subs.npc,
                                       ignore.strand = T)
e.nuc.x.gb <- GenomicRanges::intersect(e.nuc.dom.1.b, genes.subs.nuc,
                                       ignore.strand = T)

e.npc.dom.1.d <- subsetByOverlaps(nup.npc.x.elys.emb.1,
                                  c(e.npc.x.tss, e.npc.x.gb),
                                  ignore.strand = T, invert = T)
e.nuc.dom.1.d <- subsetByOverlaps(nup.nuc.x.elys.emb.1,
                                  c(e.nuc.x.tss, e.nuc.x.gb),
                                  ignore.strand = T, invert = T)

e.npc.ig.pv <- perm.test.jac.ov(e.npc.dom.1.d, intergenes)
e.nuc.ig.pv <- perm.test.jac.ov(e.nuc.dom.1.d, intergenes)

# Brain-enhancer associated genes

b.e.genes <- scan("tables/b.e.genes.id.convert.txt", what = 'character')
tss.b.e.g <- dm3.tss.gr[dm3.tss.gr$id %in% b.e.genes]
genes.b.e.g <- dm3.genes.gr[dm3.genes.gr$id %in% b.e.genes]

e.npc.b.e.g.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, tss.b.e.g)
e.npc.b.e.g.pv.ic <- perm.test.jac.ov(nup.npc.x.elys.emb.1, tss.b.e.g,
                               genome = "/home/artem/R/projects/Elys/bed/s2.inact.state.chrome.size")
e.nuc.b.e.g.pv <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, tss.b.e.g)
e.nuc.b.e.g.pv.ac <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, tss.b.e.g,
                                      genome = "/home/artem/R/projects/Elys/bed/s2.act.state.chrome.size")

subsetByOverlaps(genes.b.e.g, nup.npc.x.elys.emb.1, ignore.strand = T)

e.npc.b.e.g.2.pv <- perm.test.jac.ov(nup.npc.x.elys.emb.1, genes.b.e.g)
e.nuc.b.e.g.2.pv <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, genes.b.e.g)
e.nuc.b.e.g.2.pv.ac <- perm.test.jac.ov(nup.nuc.x.elys.emb.1, genes.b.e.g,
                                        genome = "/home/artem/R/projects/Elys/bed/s2.act.state.chrome.size")
tab <- rbind(e.npc.tss.pv, e.nuc.tss.pv,
             e.npc.gb.pv, e.nuc.gb.pv,
             e.npc.ig.pv, e.nuc.ig.pv,
             e.npc.starr.pv, e.nuc.starr.pv,
             e.npc.b.e.pv, e.nuc.b.e.pv,
             e.npc.b.e.g.pv, e.npc.b.e.g.pv.ic, e.nuc.b.e.g.pv, e.nuc.b.e.g.pv.ac,
             e.npc.b.e.g.2.pv, e.nuc.b.e.g.2.pv, e.nuc.b.e.g.2.pv.ac)

tab <- cbind(tibble(
  "A" = c("Elys x NPC*", "Elys x NUC*", "Elys x NPC*", "Elys x NUC*",
          "Elys x NPC*", "'Elys x NUC*", "Elys x NPC*", "Elys x NUC*",
          "Elys x NPC*", "Elys x NUC*", "Elys x NPC*",
          "(Elys x NPC*) x inact cols", "Elys x NUC*", "Elys x NUC* x act cols","Elys x NPC*",
          "Elys x NUC*", "Elys x NUC* x act cols"),
  "B" = c("TSS", "TSS", "Gene Bodies", "Gene Bodies",
          "Intergenic", "Intergenic", "S2 ecd enh",
          "S2 ecd enh", "Brain enh", "Brain enh", "Brain enh TSS",
          "Brain enh TSS x Inact cols", "Brain enh TSS", "Brain enh TSS x Act cols",
          "Brain enh Gene Bodies", "Brain enh Gene Bodies", "Brain enh Gene Bodies x Act cols")),
  tab)

tab$`%A`[c(1,3,5,7,9,11,12,15)] <- tab$AxB[c(1,3,5,7,9,11,12,15)]/length(nup.npc.x.elys.emb.1)*100
tab$`%A`[c(1,3,5,7,9,12,13, 15, 16)+1] <- tab$AxB[c(1,3,5,7,9,12,13, 15, 16)+1]/length(nup.nuc.x.elys.emb.1)*100
tab$`%B`[c(3,4)] <- tab$BxA[c(3,4)]/length(dm3.tss.gr)*100
write.xlsx(tab, file = "tables/elys.x.nups.overlaps.jac.xlsx", overwrite = T)

# Pie charts

elys.x.npc.plot <- data.frame(
  group = factor(c("TSS", "Gene bodies", "None"),
                 levels = c("TSS", "Gene bodies", "None")),
  value = c(tab$AxB[1]/length(nup.npc.x.elys.emb.1),
            tab$AxB[3]/length(nup.npc.x.elys.emb.1),
            1 - tab$AxB[1]/length(nup.npc.x.elys.emb.1) - tab$AxB[3]/length(nup.npc.x.elys.emb.1))
)
elys.x.nuc.plot <- data.frame(
  group = factor(c("TSS", "Gene bodies", "None"),
                 levels = c("TSS", "Gene bodies", "None")),
  value = c(tab$AxB[2]/length(nup.nuc.x.elys.emb.1),
            tab$AxB[4]/length(nup.nuc.x.elys.emb.1),
            1 - tab$AxB[2]/length(nup.nuc.x.elys.emb.1) - tab$AxB[4]/length(nup.nuc.x.elys.emb.1))
)

p1 <- ggplot(elys.x.npc.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=5) +
  scale_fill_brewer(palette="Set1", name = "Domain overlaps with")+
  ggtitle("Elys embryo x Nup98 NPC")

p2 <- ggplot(elys.x.nuc.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(x = c(1.2, 1, 1),y = ypos, label = paste0(format(prop, digits = 2),
                                                          "%")), color = "white", size=4) +
  scale_fill_brewer(palette="Set1", name = "Domain overlaps with")+
  ggtitle("Elys embryo x Nup98 NUC")

pdf("plots/elys.nup.x.genes.pie.pdf", width = 6, height = 5)
p1
p2
dev.off()


tss.elys.npc.plot <- data.frame(
  group = factor(c("TSS", "Gene body", "Nothing"),
                 levels = c("TSS", "Gene body", "Nothing")),
  value = c(tab$BxA[1]/length(dm3.tss.gr),
            tab$BxA[3]/length(dm3.tss.gr),
            1 - tab$BxA[1]/length(dm3.tss.gr) - tab$BxA[3]/length(dm3.tss.gr))
)
tss.elys.nuc.plot <- data.frame(
  group = factor(c("TSS", "Gene body", "Nothing"),
                 levels = c("TSS", "Gene body", "Nothing")),
  value = c(tab$BxA[2]/length(dm3.tss.gr),
            tab$BxA[4]/length(dm3.tss.gr),
            1 - tab$BxA[2]/length(dm3.tss.gr) - tab$BxA[4]/length(dm3.tss.gr))
)

p3 <- ggplot(tss.elys.npc.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(x = c(1, 1, 1.3),y = ypos, label = paste0(format(prop, digits = 2),
                                                          "%")), color = "white", size=3) +
  scale_fill_brewer(palette="Set1", name = "Gene overlaps\nwith domain by")+
  ggtitle("Elys embryo x Nup98 NPC")


p4 <- ggplot(tss.elys.nuc.plot %>% 
               
               arrange(desc(group)) %>%
               mutate(prop = value *100) %>%
               mutate(ypos = cumsum(prop)- 0.5*prop ),
             aes(x = "", y = prop, fill = group))+
  geom_bar(stat = "identity", width = 1, color = "white")+
  coord_polar("y", start = 0, direction = -1)+
  theme_void()+
  geom_text(aes(y = ypos, label = paste0(format(prop, digits = 2),
                                         "%")), color = "white", size=4) +
  scale_fill_brewer(palette="Set1", name = "Gene overlaps\nwith domain by")+
  ggtitle("Elys embryo x Nup98 NUC")

pdf("plots/genes.x.elys.nup.pie.pdf", width = 6, height = 5)
p3
p4
dev.off()







