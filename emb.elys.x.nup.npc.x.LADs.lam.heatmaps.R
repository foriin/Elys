library(dm3)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(genomation)

load("RData/nups.x.elys.emb.RData", verbose = T)
load("RData/lam.profiles.in.br.nrn.gl.fb.kc.RData", verbose = T)
load("RData/emb.hmm3.RData", verbose = T)
load("RData/elys.nup.npc.x.LADs.RData", verb = T)

winsor <- function(data){
  q95 <- quantile(data, 0.99)
  data[data > q95] <- q95
  data[data < -q95] <- -q95
  return(data)
}

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

elys.npc.lads <- lapply(list(npc.elys.lad.br, npc.elys.lad.fb, npc.elys.lad.glia,
                             npc.elys.lad.nrn, npc.elys.lad.kc, npc.elys.lad.emb),
                        function(nel){
                          resize(nel, width = 12000, fix = 'center')
                        })
sml.list <- lapply(1:6, function(ind){
  sml <- ScoreMatrixBin(target = profs[[ind]], windows = elys.npc.lads[[ind]],
                      strand.aware = T, bin.num = 40, weight.col = "log2damid")
  sml@.Data <- sml@.Data[order(rowSums(sml@.Data[, 19:21]),
  decreasing = F),] %>% winsor()
  sml
})
sml.sc <- scaleScoreMatrixList(ScoreMatrixList(sml.list))
sml.list <- ScoreMatrixList(sml.list, bin.num = 40)

sml1[[1]]@.Data <- sml1[[1]]@.Data[order(rowSums(sml1[[1]]@.Data[, 19:21]),
                                         decreasing = F),]
sml1[[2]]@.Data <- sml1[[2]]@.Data[order(rowSums(sml1[[2]]@.Data[, 19:21]),
                                               decreasing = F),]
sml1[[3]]@.Data <- sml1[[3]]@.Data[order(rowSums(sml1[[3]]@.Data[, 19:21]),
                                               decreasing = F),]
sml1[[4]]@.Data <- sml1[[4]]@.Data[order(rowSums(sml1[[4]]@.Data[, 19:21]),
                                               decreasing = F),]
sml1[[5]]@.Data <- sml1[[5]]@.Data[order(rowSums(sml1[[5]]@.Data[, 19:21]),
                                         decreasing = F),]
sml1[[6]]@.Data <- sml1[[6]]@.Data[order(rowSums(sml1[[6]]@.Data[, 19:21]),
                                         decreasing = F),]
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
sml2[[1]]@.Data <- sml2[[1]]@.Data[order(rowSums(sml2[[1]]@.Data[, 19:21]),
                                               decreasing = F),]
sml2[[2]]@.Data <- sml2[[2]]@.Data[order(rowSums(sml2[[2]]@.Data[, 19:21]),
                                               decreasing = F),]
sml2.sc <- scaleScoreMatrixList(sml2)
sml2.sc[[1]]@.Data <- sml2.sc[[1]]@.Data[order(rowSums(sml2.sc[[1]]@.Data[, 19:21]),
                                               decreasing = F),]
sml2.sc[[2]]@.Data <- sml2.sc[[2]]@.Data[order(rowSums(sml2.sc[[2]]@.Data[, 19:21]),
                                               decreasing = F),]

# plotMeta(sml.sc)
bluered <- colorRampPalette(c("blue3", "white", "red3"))



pdf("plots/heatmaps.lam.br.fb.gl.nrn.kc.emb.ar.elys.npc.lads.pdf", height = 10, width = 24)
par(mfrow = c(1,6))
for (i in 1:6){
  heatMatrix(sml.list[[i]], xcoords = c(-6000, 6000), col = bluered(20),
             order = F, main = names(profs)[i])
}

dev.off()

pdf("plots/heatmaps.kc.emb.ar.elys.npc.lads.pdf", height = 10, width = 8)
par(mfrow = c(1,2))
for (i in 5:6){
  heatMatrix(sml.list[[i]], xcoords = c(-6000, 6000), col = bluered(20),
             order = F, main = names(profs)[i])
}

dev.off()

pdf("plots/heatmaps.lam.kc.emb.ar.elys.npc.unsc.pdf", height = 10, width = 8)
multiHeatMatrix(sml2, xcoords = c(-6000, 6000),order = F, col = bluered(20),
                winsorize = c(0.2, 99.5), common.scale = T)
dev.off()










