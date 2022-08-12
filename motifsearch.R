# BiocManager::install("rGADEM")
# BiocManager::install(c("TFBSTools", "JASPAR2020"))
# BiocManager::install("BSgenome.Dmelanogaster.UCSC.dm3")
# BiocManager::install("seqinr")
# BiocManager::install("MotifDb")
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


load("RData/emb.lam.elys.curr.nogap.RData", verbose = T)
load("RData/nup98.damid.hmm3.nogap.RData", verbose = T)

elys.x.nup.npc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.npc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.npc <- elys.x.nup.npc[width(elys.x.nup.npc) >= 100]

elys.x.nup.nuc <- GenomicRanges::intersect(elys.emb.hmm3.ng,
                                           nup98.nuc.hmm3.uq.ng,
                                           ignore.strand = T)
elys.x.nup.nuc <- elys.x.nup.nuc[width(elys.x.nup.nuc) >= 100]
# save(elys.x.nup.npc, elys.x.nup.nuc,
#      file = "RData/elys.x.nup98.mtf.RData")

gadem.f <- function(dom, pr, score.thr = 1){
  a <- findOverlaps(dom, pr, ignore.strand = T)
  b <- lapply(unique(a@from), function(x){
    gr <- pr[a@to[a@from == x]]
    gr[which.max(gr$score)]
  })
  d <- GRangesList(b)
  e <- unlist(d)
  e <- e[e$score > score.thr]
  e <- resize(e, 1000, fix = "center")
  if (!all(grepl("^chr", seqlevels(e)))) seqlevels(e) <-
    paste0("chr", seqlevels(e))
  seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, e)
  
  GADEM(seqs, verbose = T, genome = "Dmelanogaster", pValue = 1e-8,
        minSites = length(e)/10)
}

find.jaspar <- function(gadem){
  lapply(1:length(nOccurrences(gadem)), function(ind){
    # extract the motif of interest from the GADEM object
    unknown_motif = getPWM(gadem)[[ind]]
    
    # convert the motif to a PWM matrix
    unknown_pwm   = PWMatrix(
      ID = 'unknown', 
      profileMatrix = unknown_motif
    )
    
    pwm_library = getMatrixSet(
      JASPAR2020,
      opts=list(
        collection = 'CORE',
        species    = 'Drosophila melanogaster',
        matrixtype = 'PWM'
      ))
    
    # find the most similar motif to our motif
    pwm_sim = PWMSimilarity(
      
      # JASPAR library
     test, 
      
      # out motif
      unknown_pwm,
      
      # measure for comparison
      method = 'Pearson')
    
    # extract the motif names from the pwm library
    pwm_library_list = lapply(pwm_library, function(x){
      data.frame(ID = ID(x), name = name(x))
    })
    
    # combine the list into one data frame
    pwm_library_dt = dplyr::bind_rows(pwm_library_list)
    
    # fetch the similarity of each motif to our unknown motif
    pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]
    
    # find the most similar motif in the library
    pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
    
    head(pwm_library_dt)
  })
}

find.motifDb <- function(gadem){
  motifs <- query(MotifDb, 'Dmelanogaster')
  motif.list <- motifs %>% as.list()
  pwms.q <- lapply(1:length(motif.list), function(x){
    PWMatrix("ID" = motifs@elementMetadata$providerId[x],
             profileMatrix = motif.list[[x]])
  })
  pwms.q.list <- do.call(PWMatrixList, pwms.q)
  
  
  lapply(1:length(nOccurrences(gadem)), function(ind){
    unknown_motif = getPWM(gadem)[[ind]]
    
    # convert the motif to a PWM matrix
    unknown_pwm   = PWMatrix(
      ID = 'unknown', 
      profileMatrix = unknown_motif
    )
    
    pwm_sim = PWMSimilarity(
      
      # JASPAR library
      pwms.q.list, 
      
      # out motif
      unknown_pwm,
      
      # measure for comparison
      method = 'Pearson')
    names(pwm_sim) <- ID(pwms.q.list)
    
    # extract the motif names from the pwm library
    pwm_library_list = lapply(pwms.q.list, function(x){
      data.frame(ID = ID(x), name = name(x))
    })
    
    # combine the list into one data frame
    pwm_library_dt = dplyr::bind_rows(pwm_library_list) %>% 
      cbind(., motifs@elementMetadata %>% as_tibble() %>% select(-c(1,2)))
    
    
    # fetch the similarity of each motif to our unknown motif
    pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]
    
    # find the most similar motif in the library
    pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
    
    pwm_library_dt
    
    
    
    
  })
}

gadem.npc <- gadem.f(elys.x.nup.npc, elys.emb.pr)
gadem.nuc <- gadem.f(elys.x.nup.nuc, elys.emb.pr)

save(gadem.npc, gadem.nuc, file = "RData/elys.x.nup98.gadem.RData")

nocc.npc <- nOccurrences(gadem.npc)
names(nocc.npc) <- consensus(gadem.npc)
nocc.nuc <- nOccurrences(gadem.nuc)
names(nocc.nuc) <- consensus(gadem.nuc)

motifs.npc <- find.motifDb(gadem.npc)
motifs.nuc <- find.motifDb(gadem.nuc)

names(motifs.npc) <- consensus(gadem.npc)
names(motifs.nuc) <- consensus(gadem.nuc)
write.xlsx(motifs.npc, "tables/elys.emb.x.nup.npc.gadem.motifs.search.xlsx")
write.xlsx(motifs.nuc, "tables/elys.emb.x.nup.nuc.gadem.motifs.search.xlsx")

plot(gadem.nuc[1])
plot(gadem.nuc[2])
plot(gadem.nuc[3])

pdf("plots/elys.x.nup98.npc.uq.motifs", width = 6, height = 3)
par(mfrow = c(1,3))
plot(gadem.npc[1])
plot(gadem.npc[2])
plot(gadem.npc[3])
plot(gadem.npc[4])
dev.off()

pdf("plots/elys.x.nup98.nuc.uq.motifs", width = 6, height = 3)
plot(gadem.nuc[1])
plot(gadem.nuc[2])
plot(gadem.nuc[3])
plot(gadem.nuc[4])
dev.off()

# Identify motifs

# extract the motif of interest from the GADEM object
unknown_motif = getPWM(gadem.nuc)[[3]]

# convert the motif to a PWM matrix
unknown_pwm   = PWMatrix(
  ID = 'unknown', 
  profileMatrix = unknown_motif
)

pwm_library = getMatrixSet(
  JASPAR2020,
  opts=list(
    collection = 'CORE',
    species    = 'Drosophila melanogaster',
    matrixtype = 'PWM'
  ))

# find the most similar motif to our motif
pwm_sim = PWMSimilarity(
  
  # JASPAR library
  pwm_library, 
  
  # out motif
  unknown_pwm,
  
  # measure for comparison
  method = 'Pearson')

# extract the motif names from the pwm library
pwm_library_list = lapply(pwm_library, function(x){
  data.frame(ID = ID(x), name = name(x))
})

# combine the list into one data frame
pwm_library_dt = dplyr::bind_rows(pwm_library_list)

# fetch the similarity of each motif to our unknown motif
pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]

# find the most similar motif in the library
pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]

head(pwm_library_dt)



# Random stuff
e <- elys.emb.pr[sample(1:length(elys.emb.pr), 3000)]
seqlevels(e) <- paste0("chr", seqlevels(e))
seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3, e)
names(seqs) <- paste(seqnames(e), start(e), end(e), sep = "_")
writeXStringSet(seqs, "fasta/random.300bp.seqs.fasta")


gadem.rand <- GADEM(seqs, verbose = T, genome = "Dmelanogaster", pValue = 1e-7, eValue = 10)

# S2 elys ChIP x Nup98 NUC*

summary(width(s2.chip.bed))

s2.elys.peaks.x.nup.nuc.uq <- GenomicRanges::intersect(s2.chip.bed,
                                                       nup98.nuc.hmm3.uq.ng,
                                                       ignore.strand = T)
summary(width(s2.elys.peaks.x.nup.nuc.uq))
length(s2.elys.peaks.x.nup.nuc.uq)

s2.elys.peaks.x.nup.nuc.uq <- s2.elys.peaks.x.nup.nuc.uq[width(s2.elys.peaks.x.nup.nuc.uq) > 100] %>% 
  resize(width = 100, fix = "center")
seqlevels(s2.elys.peaks.x.nup.nuc.uq) <- paste0("chr", 
                                                seqlevels(s2.elys.peaks.x.nup.nuc.uq))
seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm3,
               s2.elys.peaks.x.nup.nuc.uq)
names(seqs) <- paste(seqnames(s2.elys.peaks.x.nup.nuc.uq),
                     start(s2.elys.peaks.x.nup.nuc.uq),
                     end(s2.elys.peaks.x.nup.nuc.uq), sep = "_")
writeXStringSet(seqs, "fasta/s2.elys.x.nup.nuc.300bp.fasta")

gadem.s2 <- gadem.f(s2.elys.peaks.x.nup.nuc.uq, s2.elys.chip.bins.2)
s2.motifs <- find.motifDb(gadem.s2)


pdf("plots/s2.elys.chip.x.nup98.nuc.uq.motifs.pdf", width = 6, height = 3)
plot(gadem.s2[1])
plot(gadem.s2[2])
plot(gadem.s2[3])
plot(gadem.s2[4])
dev.off()


# Here I'll try to use library of motifs from MotifDb



motifs.s2 <- find.motifDb(gadem.s2)
write.xlsx(test, "tables/s2.elys.chip.x.nup.nuc.xlsx")

motifs.rand <- find.motifDb(gadem.rand)
lapply(motifs.rand, head)

motifs.npc <- find.motifDb(gadem.npc)
lapply(motifs.npc, head)

motifs.nuc <- find.motifDb(gadem.nuc)
lapply(motifs.nuc, head)
