library(GenomicRanges)
library(data.table)
library(dplyr)
library(tximport)
library(rtracklayer)


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
                notecol="white", notecex = 0.5, margins=c(7,7),
                main=paste(meth, "'s correlation coefficients and hierarchical clustering for\n", "'",
                           sub("\\d+\\.(.*)", "\\1", corr.desc), "'", sep=""),
                cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2),
                lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="density")
      options(warn=0)
      dev.off()
    }  
  })
}


setwd("~/IMG/Projects/Elys/kd_rnaseq/quant_euc/")


dir(pattern = "Sh")
sh1 <- fread("Sh1_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh2 <- fread("Sh2_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh3 <- fread("Sh3_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh4 <- fread("Sh4_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh5 <- fread("Sh5_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh6 <- fread("Sh6_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh7 <- fread("Sh7_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))
sh8 <- fread("Sh8_euc_quant/quant.genes.sf") %>% filter(!grepl("FBtr", Name))

s2.tpm <- data.frame(fbgn = sh1$Name,
                      dsZ1 = sh4$TPM,
                 dsE1 = sh3$TPM,
                 dsE2 = sh7$TPM)

revinf <- function(x){
  if (x == -Inf) return(Inf)
  else if (x == Inf) return(-Inf)
  else return(x)
}
s2.tpm.mol <- s2.tpm %>% filter((dsZ1 / dsE1 > 1.5 & dsZ1 / dsE2 > 1.5 & dsE1 * dsE2 * dsZ1 > 1)|
                                  (dsE1 / dsZ1 > 1.5 & dsE2 / dsZ1 > 1.5 & sum(dsE1, dsE2, dsZ1) > 0))
s2.tpm.mol <- merge(genes, s2.tpm, by.x = "id", by.y = "fbgn") %>% 
  rowwise() %>% 
  mutate(maxLogFC = ifelse(dsZ1/dsE1 > 1,
                              log2(max(revinf(dsZ1/dsE1), revinf(dsZ1/dsE2))),
                              log2(min(revinf(dsZ1/dsE1), revinf(dsZ1/dsE2))))) %>% filter(is.finite(maxLogFC), abs(maxLogFC) >= 2)
write.table(s2.tpm.mol, "kdElys.1.5.fold.change.csv", sep = ";", quote = F,
            row.names = F, dec = ",")


wilcox.test(s2.tpm.mol$dsZ1, s2.tpm.mol$dsE1)
wilcox.test(s2.tpm.mol$dsZ1, s2.tpm.mol$dsE2)
wilcox.test(s2.tpm.mol$dsE1, s2.tpm.mol$dsE2)

wilcox.test(s2.tpm$dsE1, s2.tpm$dsZ1)

elys.dom <- import.bed("/home/artem/IMG/Projects/Elys/26.03.18_elys_lam_embryos/300nt/BioHMM/ELYS.EMB.300nt.domains.bed")
lam.dom <- import.bed("/home/artem/IMG/Projects/Elys/26.03.18_elys_lam_embryos/300nt/BioHMM/LAM.EMB.300nt.domains.bed")


genes <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.full.genes.gtf") %>% 
  select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('.*gene_id "(FBgn[0-9]+)";.*', '\\1', attr), gene_name = sub('.*gene_name "([^;]+)";.*', '\\1', attr),
         tss = ifelse(strand == "+", start, end)) %>% 
  select(-attr)

s2.salmon <- merge(genes, s2.tpm, by.x = "id", by.y = "fbgn")

s2.sal.tss.gr <- GRanges(
  seqnames = Rle(s2.salmon$chr),
  ranges = IRanges(
    start = s2.salmon$tss,
    width = 1,
    names = s2.salmon$id
  )
)

write.table(s2.salmon, "kdElys.all.genes.expr.csv", sep = ";", quote = F,
            row.names = F, dec = ",")

s2.salmon <- s2.salmon %>%
  mutate(ELYS = ifelse(id %in% names(subsetByOverlaps(s2.sal.tss.gr, elys.dom)), 1, 0),
           LAM = ifelse(id %in% names(subsetByOverlaps(s2.sal.tss.gr, lam.dom)), 1, 0)) %>% 
  filter(dsZ1 > 0, dsE1 > 0, dsE2 > 0)
setwd("/home/artem/IMG/Projects/Elys/kd_rnaseq/Plots/")

pdf("elys_doms_expr.pdf", width = 12, height = 6)
  par(mfrow = c(1, 3))
  boxplot(dsZ1 ~ ELYS, s2.salmon, outline = F, main = "dsZ1", ylab = "TPM")
  boxplot(dsE1 ~ ELYS, s2.salmon, outline = F, main = "dsE1", ylab = "TPM")
  boxplot(dsE2 ~ ELYS, s2.salmon, outline = F, main = "dsE2", ylab = "TPM")
dev.off()

pdf("lam_doms_expr.pdf", width = 12, height = 6)
  par(mfrow = c(1, 3))
  boxplot(dsZ1 ~ LAM, s2.salmon, outline = F, main = "dsZ1", ylab = "TPM")
  boxplot(dsE1 ~ LAM, s2.salmon, outline = F, main = "dsE1", ylab = "TPM")
  boxplot(dsE2 ~ LAM, s2.salmon, outline = F, main = "dsE2", ylab = "TPM")
dev.off()

median(s2.salmon$dsZ1[s2.salmon$ELYS == 1])
median(s2.salmon$dsE1[s2.salmon$ELYS == 1])
median(s2.salmon$dsE2[s2.salmon$ELYS == 1])

median(s2.salmon$dsZ1[s2.salmon$ELYS == 0])
median(s2.salmon$dsE1[s2.salmon$ELYS == 0])
median(s2.salmon$dsE2[s2.salmon$ELYS == 0])

median(s2.salmon$dsZ1[s2.salmon$LAM == 1])
median(s2.salmon$dsE1[s2.salmon$LAM == 1])
median(s2.salmon$dsE2[s2.salmon$LAM == 1])

median(s2.salmon$dsZ1[s2.salmon$LAM == 0])
median(s2.salmon$dsE1[s2.salmon$LAM == 0])
median(s2.salmon$dsE2[s2.salmon$LAM == 0])
