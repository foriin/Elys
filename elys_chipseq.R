library(data.table)
library(dplyr)

setwd("~/IMG/Projects/Elys/chip_capelson/")

chip <- fread("chipseq_elys.300.bedgraph") %>% setNames(c("chr", "start", "end", "score"))
input <- fread("chipseq_input.300.bedgraph") %>% setNames(c("chr", "start", "end", "score"))

cor(chip[,4], input[, 4])

chip.norm <- chip %>% mutate(score = log2(score/input$score)) %>% filter(is.finite(score))

write.table(chip.norm, "elys_chipseq_norm.bedgraph", sep = "\t", col.names = F, row.names = F, quote = F)

chip2 <- fread("chipseq_elys.300nt.2.bedgraph") %>% setNames(c("chr", "start", "end", "score"))
input2 <- fread("chipseq_input.300nt.2.bedgraph") %>% setNames(c("chr", "start", "end", "score"))
chip.norm.2 <- chip2 %>% mutate(score = log2(score/input2$score)) %>% filter(is.finite(score))
cor(chip.norm$score, chip.norm.2$score, use = "complete.obs")

x <- merge(chip.norm, chip.norm.2, by = c("chr", "start"))

head(x)

y <- x %>% filter(score.x > 1)
cor(y$score.x, y$score.y)
library(vioplotx)

vioplotx(x$score.x, x$score.y)

write.table(chip.norm.2, "elys_chipseq_norm.2.bedgraph", sep = "\t", col.names = F, row.names = F, quote = F)


chip.sum <- chip %>% mutate(score = chip$score + chip2$score)

input.sum <- input %>% mutate(score = input$score + input2$score)

chip.sum.norm <- chip.sum %>% mutate(score = log2(score/input.sum$score)) %>% filter(is.finite(score))

write.table(chip.sum.norm, "elys_chipseq_norm.sum.bedgraph", sep = "\t", col.names = F, row.names = F, quote = F)

