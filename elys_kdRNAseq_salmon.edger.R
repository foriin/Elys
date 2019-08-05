library(GenomicRanges)
library(data.table)
library(dplyr)
library(stringr)
library(gplots)
library(tximport)
library(edgeR)

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

emb.tpm <- data.frame(fbgn = sh1$Name,
                      dsZ1 = sh4$TPM,
                 dsZ2 = sh8$TPM,
                 dsE1 = sh3$TPM,
                 dsE2 = sh7$TPM)

emb.sh <- c("Sh4_euc_quant", "Sh8_euc_quant", "Sh3_euc_quant", "Sh7_euc_quant")
emb.files <- file.path(emb.sh, "quant.genes.sf")
names(emb.files) <- c("dsZ1", "dsZ2", "dsE1", "dsE2")
txi.emb <- tximport(emb.files, "salmon", txOut = T)

head(txi.emb$counts)

dgList <- DGEList(counts = txi.emb$counts, genes = row.names(txi.emb$counts))
dgList
head(dgList$counts)

dg.cpm <- cpm(dgList)
summary(dg.cpm)
# Filtering
ccheck <- dg.cpm > 1
summary(ccheck)

keep <- which(rowSums(ccheck) >= 2)
dgList <- dgList[keep,]
summary(cpm(dgList))
# Normalise 
dgList <- calcNormFactors(dgList, method = "TMM")

plotMDS(dgList, dim.plot = c(1,3))

anno <- fread("~/IMG/data/dmel/gene_list/Drosophila_melanogaster.BDGP5.78.gtf") %>% filter(V3 == "gene")

fbgn.cg <- str_split_fixed(anno$V9, ";", n = 5) %>% as.data.frame() %>% select(1,3,5) 
x <- apply(fbgn.cg, 2, function(x) gsub('.*"([^"]+)".*', '\\1', x))
fbgn.cg <- as.data.frame(x) %>% setNames(c("FBgn", "gene_name", "type"))

shall <- merge(fbgn.cg, shall, by = "FBgn")

MainCorrelations(shall[, 4:12], corr.desc = "TPMs", file.p = ".")

setwd("~/IMG/Projects/Elys/kd_rnaseq/quant_norib/")

files <- file.path(dir(pattern = "norib"), "quant.genes.sf")
names(files) <- paste0("Sh", 1:8)
txi.salmon <- tximport(files[c(1,2,5,6)], type = "salmon", txOut = T)
head(txi.salmon$counts)

group <- factor(c(1,1,2,2))
y <- DGEList(counts=txi.salmon$counts, group = group)

y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)
topTags(qlf)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
