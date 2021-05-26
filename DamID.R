library(dm3)
library(GenomicAlignments)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(gplots)

MakeCorMat <- function(data, cormethod, use.vals = "complete.obs"){
  cormat <- as.matrix(data) %>%
    cor(method = cormethod, use = use.vals) %>% 
    round(digits = 2)
}

MainCorrelations <- function(data, use.opt="everything", corr.desc, createPDF=T, file.p, ...) {  
  lapply(c("spearman", "pearson"), function(meth){
    cors <- MakeCorMat(data, meth)
    assign(paste0(corr.desc, ".", meth, ".cor"), cors)
    if (createPDF == T){
      options(warn=-1)
      pdf(file=file.path(file.p, paste0(corr.desc,  "_", meth, "_correlation_heatmap", ".pdf")), width=12, height=12)
      heatmap.2(cors, col=greenred(200), breaks=seq(from=-1, to=1, by=0.01),
                Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=cors,
                notecol="white", notecex = 1.5, margins=c(7,7),
                main=paste(meth, "'s correlation coefficients and hierarchical clustering for\n", "'",
                           sub("\\d+\\.(.*)", "\\1", corr.desc), "'", sep=""),
                cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2),
                lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="density")
      options(warn=0)
      dev.off()
    }  
  })
}

filt.dam <- dm3.genes.gr[grepl("Hsp70|^w$",dm3.genes.gr$gene_name)]
filt.lam <- dm3.genes.gr[grepl("Hsp70|Lam$|^w$",dm3.genes.gr$gene_name)]
filt.elys <- dm3.genes.gr[grepl("Hsp70|CG14215|^w$",dm3.genes.gr$gene_name)]

filt.l <- list("D" = filt.dam,
               "DL" = filt.lam,
               "DE" = filt.elys)
filt.all <- dm3.genes.gr[grepl("Hsp70|CG14215|Lam$|^w$",dm3.genes.gr$gene_name)]


covs <- new.env()
bam.dir <- "/home/artem/IMG/Projects/Elys/DamID_25.12.20/dm3.map"
# test <- readGAlignments("/home/artem/IMG/Projects/Elys/DamID_25.12.20/dm3.map/IED1.dm3.map.bam")
dir(bam.dir, pattern = "bam$")

for (i in dir(bam.dir, pattern = "bam$")){
  bam <- readGAlignments(file.path(bam.dir, i)) 
  # filt.ind <- grep(sub(".{2}(D[LE]?)[12].*", "\\1", i), names(filt.l))[1]
  # bam <- subsetByOverlaps(GRanges(bam), filt.l[[filt.ind]], maxgap = 510,
                          # invert = T, ignore.strand = T)
  bam <- subsetByOverlaps(GRanges(bam[njunc(bam) == 0]), filt.all, maxgap = 510,
                          invert = T, ignore.strand = T)
  bins <- tileGenome(seqinfo(bam), tilewidth = 500, cut.last.tile.in.chrom = T)
  assign(sub(".dm3.map.bam", ".rpm", i),
         binnedAverage(bins, coverage(bam)/length(bam)*1e6, varname = "score"),
         envir = covs)
}

score.mat <- sapply(ls(envir = covs), function(x){
  get(x, envir = covs)$score
})

score.mat.euc <- sapply(ls(envir = covs), function(x){
  subsetByOverlaps(get(x, envir = covs),
                   euc.gr)$score
})

MainCorrelations(score.mat, corr.desc = "DamID.RPM", file.p = "plots")
MainCorrelations(score.mat.euc, corr.desc = "DamID.RPM.euc", file.p = "plots")
dir.create("plots/repcor")
for (n in seq(1, by = 2, length.out = 7)){
  namens <- sub("1.rpm", "", colnames(score.mat)[n])
  png(paste0("plots/repcor/", namens, ".repcor.png"))
  plot(log2(score.mat[,n]), log2(score.mat[,n+1]), type = "p",
       main = namens)
  text(-10, 5, paste(
    "Corr Sp =", format(cor(score.mat[,n], score.mat[,n+1]), digits = 2)))
  abline(coef = c(0,1), col = "red")
  dev.off()
  
}


for (cover in ls(env = covs)){
  export.bedGraph(
    subsetByOverlaps(
      get(cover, env = covs),
      euc.gr
    ), con = paste0("bed/", cover, ".bedgraph")
  )
}

test <- df.from.GRanges(covs$IEDL1.rpm)
# test$score <- covs$IEDL1.rpm$score
test <- cbind(test, score.mat)

covs.n <- lapply(ls(env = covs, pattern = "D[LE][12]"), function(rpm){
  covv <- get(rpm, envir = covs)
  covv$score <- log2(covv$score/(get(
    sub("(..)[LE]([12].rpm)", "\\1\\2", rpm),
    env = covs
  )$score))
  covv[is.finite(covv$score)]
}
)

names(covs.n) <- sub(".rpm", ".norm", ls(env = covs, pattern = "D[LE][12]"))

for (norm in names(covs.n)){
  export.bedGraph(
    subsetByOverlaps(
      covs.n[[norm]],
      euc.gr
    ), con = paste0("bed/", norm, ".bedgraph")
  )
}

