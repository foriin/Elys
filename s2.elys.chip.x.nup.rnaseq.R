library(dm3)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(openxlsx)

load("RData/kc.nup98.domains.RData", verbose = T)
load("RData/s2.elys.chip.RData", verbose = T)
load("RData/s2.tpms.RData", verbose = T)
load("RData/s2.lam.elys.deseq.RData", verbose = T)

# Subset only those ELYS peaks that overlap with either NPC fraction of 
# Nup198 or nucleoplasmic fraction
elys.x.nup.npc <- subsetByOverlaps(s2.chip.bed, kc.nup.npc, ignore.strand = T)
elys.x.nup.nuc <- subsetByOverlaps(s2.chip.bed, kc.nup.nuc, ignore.strand = T)

# Find genes that overlap with ranges found above by their TSSs
genes.x.e.n.npc <- subsetByOverlaps(dm3.tss.gr,
                                    elys.x.nup.npc, ignore.strand = T)$id
genes.x.e.n.nuc <- subsetByOverlaps(dm3.tss.gr,
                                    elys.x.nup.nuc, ignore.strand = T)$id
id.1 <- setdiff(genes.x.e.n.nuc, genes.x.e.n.npc)
id.2 <- setdiff(genes.x.e.n.npc, genes.x.e.n.nuc)
# See their comparative expression levels
boxplot(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.1],
        s2.elys.TPM$ElysTPMav[s2.elys.TPM$id %in% id.1],
        outline = F)
wilcox.test(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.1],
            s2.elys.TPM$ElysTPMav[s2.elys.TPM$id %in% id.1],
            alt = "l")
wilcox.test(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.2],
            s2.elys.TPM$ElysTPMav[s2.elys.TPM$id %in% id.2],
            alt = "l")

# Check genes that are on X and autosomes separately

id.1.x <- base::intersect(id.1,
                  dm3.genes %>% filter(chr == "X") %>% pull(id))
id.1.a <- base::intersect(id.1,
                  dm3.genes %>% filter(chr %in% c("2L", "2R", "3L", "3R")) %>%
                    pull(id))

boxplot(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.1.a],
        s2.elys.TPM$ElysTPMav[s2.elys.TPM$id %in% id.1.a],
        outline = F, main = "Elys x Nup98 nuc on A",
        ylab = "TPM", names = c("WT", "Elys KD"))
boxplot(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.1.x],
        s2.elys.TPM$ElysTPMav[s2.elys.TPM$id %in% id.1.x],
        outline = F, main = "Elys x Nup98 nuc on X",
        ylab = "TPM", names = c("WT", "Elys KD"))

# Split genes by their expression level
s2.elys.tpm.id.1 <- s2.elys.TPM %>% filter(id %in% id.1) %>% 
  mutate(EKD_norm_WT = ElysTPMav/ZTPMav,
         TPM_cut = cut(ZTPMav, c(0,1,10,100, max(ZTPMav)), include.lowest = T,
         labels = c("1", "2", "3", "4")))

boxplot(s2.elys.tpm.id.1$EKD_norm_WT, outline = F)
wilcox.test(s2.elys.tpm.id.1$EKD_norm_WT, mu = 1, alt = "g")


tpm.cut.E.Z.pv <- sapply(c("1", "2", "3", "4"), function(lev){
  format(
    wilcox.test(
      s2.elys.tpm.id.1 %>% filter(TPM_cut == lev) %>% 
        pull(ZTPMav),
      s2.elys.tpm.id.1 %>% filter(TPM_cut == lev) %>% 
        pull(ElysTPMav)
    )$p.value, digits = 3
  ) 
})

s2.elys.tpm.id.2 <- s2.elys.tpm.id.1 %>% mutate(EKD_norm_WT = ElysTPMav/ZTPMav) %>% 
  group_by(TPM_cut) %>% 
  summarize(ElysKD_TPMav = mean(ElysTPMav),
            WT_TPMav = mean(ZTPMav),
            EKD_norm_WTav = mean(EKD_norm_WT, na.rm = T, trim = 0.01),
            n = n()) %>% 
  mutate(pv = tpm.cut.E.Z.pv,
         TPM_cut = c("0-1", "1-10", "10-100", ">100"))

boxplot(EKD_norm_WT ~ TPM_cut, data = s2.elys.tpm.id.1, outline = F)
tpm.norm.cut.E.Z.pv <- sapply(c("1", "2", "3", "4"), function(lev){
  format(
    wilcox.test(
      s2.elys.tpm.id.1 %>% filter(TPM_cut == lev) %>% 
        pull(EKD_norm_WT), mu = 1, alt = "g",
    )$p.value, digits = 3
  ) 
})
s2.elys.tpm.id.3 <- s2.elys.tpm.id.1 %>% 
  group_by(TPM_cut) %>% 
  summarize(ElysKD_TPMav = mean(ElysTPMav),
            WT_TPMav = mean(ZTPMav),
            EKD_norm_WTav = mean(EKD_norm_WT, na.rm = T, trim = 0.01),
            n = n()) %>% 
  mutate(pv = tpm.norm.cut.E.Z.pv,
         TPM_cut = c("0-1", "1-10", "10-100", ">100"))

# Ranges of genes, differentially expressed upon Elys KD
upreg <- res.df.elyskd %>% filter(padj < 0.05, log2FoldChange > log2(1.5))
downreg <- res.df.elyskd %>% filter(padj < 0.05, log2FoldChange < log2(1/1.5))

up.gr <- dm3.genes.gr[dm3.genes.gr$id %in% upreg$id]
down.gr <- dm3.genes.gr[dm3.genes.gr$id %in% downreg$id]

# See if these genes are juxtaposed with subsets of chip peaks found above:

id.1 <- subsetByOverlaps(up.gr, elys.x.nup.npc, ignore.strand = T)$id
id.2 <- subsetByOverlaps(up.gr, elys.x.nup.nuc, ignore.strand = T)$id

id.3 <- subsetByOverlaps(down.gr, elys.x.nup.npc, ignore.strand = T)$id
id.4 <- subsetByOverlaps(down.gr, elys.x.nup.nuc, ignore.strand = T)$id

boxplot(s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.1],
        s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.2],
        s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.3],
        s2.elys.TPM$ZTPMav[s2.elys.TPM$id %in% id.4],
        outline = F)
plot(s2.elys.TPM$ZTPMav, s2.elys.TPM$ElysTPMav)
abline(a = 1, b = 1, col = 'red')


elys.x.nup.nuc.w <- sum(width(GenomicRanges::intersect(s2.chip.bed, 
                                                       kc.nup.nuc,
                                                       ignore.strand = T)))

elys.x.nup.nuc.pv <- perm.test.length(s2.chip.bed,
                       kc.nup.nuc, contr.val = elys.x.nup.nuc.w)

elys.x.nup.npc.w <- sum(width(GenomicRanges::intersect(s2.chip.bed, 
                                                       kc.nup.npc,
                                                       ignore.strand = T)))
elys.x.nup.npc.pv <- perm.test.length(s2.chip.bed,
                                      kc.nup.npc, contr.val = elys.x.nup.npc.w)
save(elys.x.nup.nuc, elys.x.nup.npc,
     genes.x.e.n.npc, genes.x.e.n.nuc,
     id.1, id.1.a, id.1.x, id.2,
     s2.elys.tpm.id.1,
     s2.elys.tpm.id.2,
     s2.elys.tpm.id.3,
     file = "RData/s2.elys.x.nup.RData")

# 