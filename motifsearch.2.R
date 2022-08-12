library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(seqinr)
library(rGADEM)
library(BSgenome.Dmelanogaster.UCSC.dm3)
library(JASPAR2020)
library(TFBSTools)
library(MotifDb)
library(openxlsx)
library(gplots)

load("RData/elys.x.nup98.gadem.RData", verbose = T)

comp.gadem <- function(gadem.1, gadem.2){
  lib <- lapply(1:length(nOccurrences(gadem.1)), function(x){
    PWMatrix(
      ID = names(getPWM(gadem.1))[x],
      profileMatrix = getPWM(gadem.1[[x]])
    )
  })%>% do.call(PWMatrixList, .)
  
  qlib <- lapply(1:length(nOccurrences(gadem.2)), function(x){
    PWMatrix(
      ID = names(getPWM(gadem.2))[x],
      profileMatrix = getPWM(gadem.2[[x]])
    )
  })
  
  ret <- sapply(qlib, function(set){
    sim <- PWMSimilarity(lib, set, method = "Euclidean")
    names(sim) <- names(getPWM(gadem.1))
    sim 
  }) %>% as.data.frame() %>% setNames(names(getPWM(gadem.2)))
  ret
}


npc.vs.nuc <- comp.gadem(gadem.npc, gadem.nuc)
nuc.vs.npc <- comp.gadem(gadem.nuc, gadem.npc)
pdf("plots/motif.heatmap.pdf" )
heatmap.2(as.matrix(npc.vs.nuc), Rowv = F, Colv = F,
          dendrogram = "none", trace = "none", margins = c(10,12), 
          cexCol = 0.7, cexRow = 0.7, key = T, col = heat.colors(n = 10, alpha = 1, rev = T))
dev.off()

heatmap(as.matrix(npc.vs.nuc), Rowv = F, Colv = F, margins = c(7,7) )





test <- lapply(1:length(nOccurrences(gadem.npc)), function(x){
  PWMatrix(
    ID = names(getPWM(gadem.npc))[x],
    profileMatrix = getPWM(gadem.npc[[x]])
  )
}) %>% do.call(PWMatrixList, .)
test <- do.call(PWMatrixList,test)
tset <- lapply(1:length(nOccurrences(gadem.nuc)), function(x){
  PWMatrix(
    ID = names(getPWM(gadem.nuc))[x],
    profileMatrix = getPWM(gadem.nuc[[x]])
  )
})

sim.1 <- sapply(tset, function(set){
  sim <- PWMSimilarity(test, set, method = "Euclidean")
  names(sim) <- names(getPWM(gadem.npc))
  sim
})
