library(dplyr)
library(ggplot2)

load("RData/tximport.RData")

# Merge difseq data and TPM data
ab <- txi.lartest$abundance %>% as.data.frame()
tpms <- data.frame(
  id = row.names(ab),
  WT = rowMeans(ab[, 1:3]),
  ElysKD = rowMeans(ab[, 4:6]),
  LamKD = rowMeans(ab[, 7:9])
)

tpms <- merge(genes, tpms, by = "id") %>% 
  filter(chr %in% c("chr2L", "chr2R", "chr3L", "chr3R",
                    "chrX")) %>% 
  mutate(chrom = ifelse(chr == "chrX", "X", "A"),
         chrom = factor(chrom, levels = c("X", "A")))
tpms[, 8:10][tpms[, 8:10] == 0] <- min(tpms[, 8:10][tpms[, 8:10] > 0])
tpms.m <- melt(tpms, measure.vars = c("WT", "ElysKD", "LamKD"))
tpms.chrom.sp <- split(tpms.m, list(tpms.m$variable, tpms.m$chrom))
# Wilcox pairwise tests
mapply(function(x,y) format(wilcox.test(tpms.chrom.sp[[x]]$value,
                                        tpms.chrom.sp[[y]]$value)$p.value, digits = 3),
       1:6,
       c(2,3,1, 5, 6, 4))

# Now for the ratios of expressions

# Ratios of gene expressions


tpms.1 <- tpms %>% mutate(el.wt = ElysKD/WT,
                          lam.wt = LamKD/WT)

boxplot(lam.wt ~ chrom, data = tpms.1,
        outline = F)

boxplot(el.wt ~ chrom, data = tpms.1,
        outline = F)
tapply(tpms.1$el.wt, tpms.1$chrom, median, na.rm = T)

boxplot(lam.wt ~ chrom, data = tpms.1,
        outline = F)
tapply(tpms.1$lam.wt, tpms.1$chrom, median, na.rm = T)

wilcox.test(el.wt ~ chrom, data = tpms.1, alt = "l")
wilcox.test(lam.wt ~ chrom, data = tpms.1, alt = "l")
boxplot(lam.wt ~ chrom, data = tpms.1,
        outline = F)
wilcox.test(lam.wt ~ chrom, data = tpms.1)
wilcox.test(el.wt ~ chrom, data = tpms.1)

# Do the same only for genes that are in any experiment have more than 1 TPM

tpms.3 <- tpms.1 %>% filter(WT > 1, ElysKD > 1, LamKD > 1)
tpms.3.m <- melt(tpms.3, measure.vars = c("WT", "ElysKD", "LamKD"))
tpms.3.chrom.sp <- split(tpms.3.m, list(tpms.3.m$variable, tpms.3.m$chrom))




# Subset genes in elys KD that are differentially expressed via padj < 0.05
# and their TPM ratio between elys KD and WT is more than 1.5 or less than 1/1.5
res.df.elys.sig <- res.df.elys %>% filter(padj < 0.05)
tpms.el.s <- tpms %>% filter(abs(log2(ElysKD/WT)) > 1,
                             id %in% res.df.elys.sig$id) %>% 
  mutate(el.wt = ElysKD/WT)
tpms.el.s.m <- melt(tpms.el.s, measure.vars = c("WT", "ElysKD"))
tapply(tpms.el.s$el.wt, tpms.el.s$chrom, median)

# Do the same for lamins KD
res.df.lam.sig <- res.df.lam %>% filter(padj < 0.05)
tpms.lam.s <- tpms %>% filter(abs(log2(LamKD/WT)) > 1,
                             id %in% res.df.lam.sig$id) %>% 
  mutate(lam.wt = LamKD/WT)
tpms.lam.s.m <- melt(tpms.lam.s, measure.vars = c("WT", "LamKD"))


# Dot plots

# Elys KD


pdf("plots/difex.elys.scatter.pdf")
ggplot(tpms %>% filter(id %in% res.df.elys.sig$id,
                       WT > 0, ElysKD > 0) %>%
         mutate(chrom = ifelse(chr == "chrX", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(ElysKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()
dev.off()

# Ugly

# Lam KD
pdf("plots/difex.lam.scatter.pdf")
ggplot(tpms %>% filter(id %in% res.df.lam.sig$id, LamKD > 0, WT >0) %>%
         mutate(chrom = ifelse(chr == "chrX", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(LamKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()
dev.off()


# Check SpC-specific genes


spc.genes <- scan("/home/artem/IMG/Projects/LAM.SpG.SpC/RNA-seq/sclspec.fbgns.2.txt",
                  what = "char")

tpms %>% filter(id %in% spc.genes)

ggplot(tpms %>% filter(id %in% spc.genes, ElysKD > 0, WT >0) %>%
         mutate(chrom = ifelse(chr == "chrX", "X", "A")) %>% 
         arrange(chrom),
       aes(x = log2(WT), y = log2(ElysKD), col = chrom))+
  geom_point()+
  geom_abline()+
  theme_bw()

boxplot((tpms %>% filter(id %in% spc.genes))$ElysKD,
        (tpms %>% filter(id %in% spc.genes))$WT,
        outline = F)
table((tpms %>% filter(id %in% spc.genes))$chr)

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% c("chr2L", "chr2R", "chr3L", "chr3R")))$WT,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "chrX"))$WT,
        outline = F)

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% c("chr2L", "chr2R", "chr3L", "chr3R")))$ElysKD,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "chrX"))$ElysKD,
        outline = F)

boxplot((tpms %>% filter(id %in% spc.genes,
                         chr %in% c("chr2L", "chr2R", "chr3L", "chr3R")))$LamKD,
        (tpms %>% filter(id %in% spc.genes,
                         chr == "chrX"))$LamKD,
        outline = F)

# Check ratios for spc specific genes

tpms.2 <- tpms %>% filter(id %in% spc.genes) %>% 
  mutate(el.wt = ElysKD/WT,
                          lam.wt = LamKD/WT,
                          chrom = ifelse(chr == "chrX", "X", "A"),
         chrom = factor(chrom, levels = c("X", "A")))
tpms.2.m <- melt(tpms.2, measure.vars = c("WT", "ElysKD", "LamKD"))

tpms.2.sp <- split(tpms.2.m, list(tpms.2.m$variable, tpms.2.m$chrom))
mapply(function(x,y) format(wilcox.test(x, y)$p.value, digits = 2),
       lapply(tpms.chrom.sp[c(1,3,5)], function(x) x$value),
       lapply(tpms.chrom.sp[c(2,4,6)], function(x) x$value))
tapply(tpms.2$el.wt, tpms.2$chrom, median, na.rm = T)

boxplot(lam.wt ~ chrom, data = tpms.2,
        outline = F)
tapply(tpms.2$lam.wt, tpms.2$chrom, median, na.rm = T)

wilcox.test(lam.wt ~ chrom, data = tpms.2, alt = "l")
wilcox.test(el.wt ~ chrom, data = tpms.2, alt = "l")


# Save data
save(tpms, tpms.m, spc.genes,
     tpm.sig,tpms.chrom.sp,
     tpms.el.s, tpms.el.s.m,
     tpms.lam.s, tpms.lam.s.m,
     res.df.elys.sig, res.df.lam.sig,
     tpms.1, tpms.2, tpms.2.m, tpms.2.sp,
     tpms.3, tpms.3.m, tpms.3.chrom.sp,
     file = "RData/tpms.analysis.RData")
