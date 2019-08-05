library(GenomicRanges)
library(data.table)
library(dplyr)

rm(list=ls())
setwd("~/IMG/Projects/Elys/26.03.18_elys_lam_embryos/")

elys.hmm.300 <- fread("./300nt/BioHMM/ELYS.EMB.300nt.domains.bed", skip = 1)
lam.hmm.300 <- fread("./300nt/BioHMM/LAM.EMB.300nt.domains.bed", skip = 1)

elys.300.gr <- GRanges(
  seqnames = Rle(elys.hmm$V1),
  ranges = IRanges(
    start = elys.hmm$V2,
    end = elys.hmm$V3
  )
)

lam.300.gr <- GRanges(
  seqnames = Rle(lam.hmm$V1),
  ranges = IRanges(
    start = lam.hmm$V2,
    end = lam.hmm$V3
  )
)

elys.w.300 <- sum(width(elys.gr))
lam.w.300 <- sum(width(lam.gr))

lam.x.elys.300 <- GenomicRanges::intersect(elys.gr, lam.gr)
elys.x.lam.w.300 <- sum(width(lam.x.elys))

elys.x.lam.w.300/lam.w.300
elys.x.lam.w.300/elys.w.300


# Shuffle regions

bedTools.shuffle.jac <- function(bed.1, bed.2, shuf = F, opt.string="-chrom"){
  
  bed.file.1 <- tempfile()
  bed.file.2 <- tempfile()
  
  shuf.1 <- tempfile()
  shuf.2 <- tempfile()
  
  jac <- tempfile()
  
  options(scipen = 99)
  
  write.table(bed.1, file = bed.file.1, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(bed.2, file = bed.file.2, quote = F, sep = "\t", col.names = F, row.names = F)
  if (shuf){
    command = paste("bedtools shuffle -i", bed.file.1,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.1, ";",
                    "bedtools shuffle -i", bed.file.2,
                    "-g /home/artem/IMG/data/dmel/Genome/dm3.genome", opt.string, "|",
                    "bedtools sort -i - >", shuf.2, ";",
                    "bedtools jaccard -a", shuf.1, "-b", shuf.2, ">", jac)
    # cat(command, "\n")
  } else {
    command = paste("bedtools jaccard -a", bed.file.1, "-b", bed.file.2, ">", jac)
  }
  try(system(command), silent = T)
  
  res=read.table(jac, header = T)
  unlink(bed.file.1); unlink(bed.file.2); unlink(shuf.2); unlink(shuf.1); unlink(jac)
  return(res$jaccard)
}


# load("jac.perm.test.hp1.lam.RData")

chr.len <- c("chr2L" = 23011544,
             "chr2R" = 21146708,
             "chr3L" = 24543557,
             "chr3R" = 27905053,
             "chr4" = 1351857,
             "chrX" = 22422827,
             "chr2LHet" = 368872,
             "chr2RHet" = 3288761,
             "chr3LHet" = 2555491,
             "chr3RHet" = 2517507,
             "chrXHet" = 204112,
             "chrYHet" = 347038)

 euc <- c("2L", "2R", "3L", "3R", "X")
 
 x <- tryCatch(bedTools.shuffle.jac(elys.hmm, elys.hmm, shuf = T),
               error = function(e) NA)

jac.300 <- bedTools.shuffle.jac(elys.hmm.300, lam.hmm.300)
print(jac)
jac.shuf.300 <- sapply(1:12000, function(i){
  tryCatch(bedTools.shuffle.jac(elys.hmm.300, elys.hmm.300, shuf = T),
           error = function(e) NA)
})
sum(!is.na(jac.shuf.300))
pval <- sum(jac.shuf >= jac, na.rm = T)/10000

print(pval)


elys.hmm.500 <- fread("500nt/BioHMM/ELYS.EMB.500nt.domains.bed", skip = 1)
lam.hmm.500 <- fread("500nt/BioHMM/LAM.EMB.500nt.domains.bed", skip = 1)

jac.500 <- bedTools.shuffle.jac(elys.hmm.500, lam.hmm.500)

jac.shuf.500 <- sapply(1:12000, function(i){
  tryCatch(bedTools.shuffle.jac(elys.hmm.500, elys.hmm.500, shuf = T),
           error = function(e) NA)
})

sum(!is.na(jac.shuf.500))
pval <- sum(jac.shuf.500 >= jac.500, na.rm = T)/10000

median(jac.shuf.500, na.rm = T)
median(jac.shuf, na.rm = T)
