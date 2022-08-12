library(dm3)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(genomation)

load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/lam.profiles.in.br.nrn.gl.fb.kc.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)

profs[[6]] <- lam.emb.pr
profs[[6]]@metadata <- list(
  protein = "Lam",
  title = "Lam DamID in embryos",
  y_lab = "log2(Dam/Dam-Lam)"
)

profs <- lapply(profs, function(pr){
  pr$log2damid <- pr$score
  pr[!is.na(pr$log2damid)]
})
names(profs) <- c("Whole brains", "Fat bodies", "Glia", "Neurons", "Kc167", "Embryos")

elys.npc.center <- resize(nup.npc.x.elys.emb.1, width = 12000, fix = 'center')


sml1 <- ScoreMatrixList(targets = profs[1:4], windows = elys.npc.center,
                       strand.aware = T,weight.col = "log2damid", bin.num = 40)
sml1.sc <- scaleScoreMatrixList(sml1)
sml1.sc[[1]]@.Data <- sml1.sc[[1]]@.Data[order(rowSums(sml1.sc[[1]]@.Data[, 19:21]),
                                             decreasing = F),]
sml1.sc[[2]]@.Data <- sml1.sc[[2]]@.Data[order(rowSums(sml1.sc[[2]]@.Data[, 19:21]),
                                             decreasing = F),]
sml1.sc[[3]]@.Data <- sml1.sc[[3]]@.Data[order(rowSums(sml1.sc[[3]]@.Data[, 19:21]),
                                             decreasing = F),]
sml1.sc[[4]]@.Data <- sml1.sc[[4]]@.Data[order(rowSums(sml1.sc[[4]]@.Data[, 19:21]),
                                             decreasing = F),]

sml2 <- ScoreMatrixList(targets = profs[5:6], windows = elys.npc.center,
                        strand.aware = T,weight.col = "log2damid", bin.num = 40)
sml2.sc <- scaleScoreMatrixList(sml2)
sml2.sc[[1]]@.Data <- sml2.sc[[1]]@.Data[order(rowSums(sml2.sc[[1]]@.Data[, 19:21]),
                                               decreasing = F),]
sml2.sc[[2]]@.Data <- sml2.sc[[2]]@.Data[order(rowSums(sml2.sc[[2]]@.Data[, 19:21]),
                                               decreasing = F),]

# plotMeta(sml.sc)
bluered <- colorRampPalette(c("blue3", "white", "red3"))



pdf("plots/heatmaps.lam.br.fb.gl.nrn.ar.elys.npc.pdf", height = 10, width = 16)
multiHeatMatrix(sml1.sc, xcoords = c(-6000, 6000),order = F, col = bluered(20),
                winsorize = c(1, 99), common.scale = T)
dev.off()

pdf("plots/heatmaps.lam.kc.emb.ar.elys.npc.pdf", height = 10, width = 8)
multiHeatMatrix(sml2.sc, xcoords = c(-6000, 6000),order = F, col = bluered(20),
                winsorize = c(1, 99), common.scale = T)
dev.off()










