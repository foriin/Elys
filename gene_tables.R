library(dplyr)
library(data.table)
library(GenomicRanges)
library(openxlsx)
library(tibble)

# Load genes (dm3 r5.57)

genes <- fread("~/Database/dmel/Annotations/dmel-gene-r5.57.gff") %>% 
  dplyr::select(c(1, 4, 5, 7, 9)) %>% 
  setNames(c("chr", "start", "end", "strand", "attr")) %>% 
  mutate(id = sub('^ID=(FBgn[0-9]+);.*', '\\1', attr),
         chr = paste0("chr", chr), gene_name = sub('.*Name=([^;]+);.*', '\\1', attr),
         tss = ifelse(strand == "+", start, end)) %>% 
  dplyr::select(-attr)

# Load testis- and bam-specific genes, define spermatocyte-spec genes
test.spec <- scan("~/IMG/Projects/LAM.SpG.SpC/RNA-seq/testspec.fbgns.3.txt",
                  what = character(), sep = '\n')


ubiq.genes <- scan("~/IMG/Projects/LAM.SpG.SpC/RNA-seq/ubiq.chiant.genes.upd.txt",
                   what = character(), sep = '\n')
# Load Laktionov RNA-seq data in bam and wt testes

load("/home/artem/R/projects/Lam_sperm/tximport.RData")

bam.tpms <- data.frame(
  id = rownames(txi.bam.wt$abundance),
  TPM_bam = rowMeans(txi.bam.wt$abundance[, 1:2])
) %>% remove_rownames()

bam.tpms <- merge(genes, bam.tpms, by = "id")

# Load our RNA-seq in spc

load("RData/tximport.RData")

lartest.tpms <- data.frame(
  id = rownames(txi.lartest$abundance),
  TPM_SpC = rowMeans(txi.lartest$abundance[, 1:3],),
  TPM_SpC_ElysKD = rowMeans(txi.lartest$abundance[, 4:6],),
  TPM_SpC_LambcKD = rowMeans(txi.lartest$abundance[, 7:9],)
) %>% remove_rownames()

tpms.all <- merge(bam.tpms, lartest.tpms, by = "id") %>%
  mutate(testspec = ifelse(id %in% test.spec, 1, 0),
         bamexp = ifelse(id %in% test.spec & TPM_bam > 1, 1, 0),
         spcspec = ifelse(id %in% test.spec & TPM_bam < 1 & TPM_SpC > 1, 1, 0),
         ubiq = ifelse(id %in% ubiq.genes & TPM_bam > 1 & TPM_SpC > 1, 1, 0))

# tpms.f <- tpms.all %>% filter(testspec == 1, !(bamspec < 1 & spcspec < 1) )
# test.spec <- tpms.f$id
# cat(test.spec, file = "~/IMG/Projects/LAM.SpG.SpC/RNA-seq/testspec.fbgns.3.txt",
#     sep = "\n")


# Load difex data

load("RData/lartest.DESeq.RData")

tpms.all <- tpms.all %>% 
  mutate(ElysKD_difex = ifelse(abs(log2(TPM_SpC_ElysKD/TPM_SpC)) > 1 &
                                 id %in% res.df.elys.sig$id, 1, 0),
         LambcKD_difex = ifelse(abs(log2(TPM_SpC_LambcKD/TPM_SpC)) > 1 &
                                 id %in% res.df.lam.sig$id, 1, 0))



write.xlsx(tpms.all, file = "./tables/SpG_and_SpC_tpms_and_spec_genes_subsets.xlsx")
