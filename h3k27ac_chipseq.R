library(data.table)
library(dplyr)

setwd("~/IMG/Projects/Elys/chip_capelson/")

h3k27.1 <- fread("chipseq_h3k27ac.300nt.1.bedgraph") %>% setNames(c("chr", "start", "end", "score"))
h3k27.2 <- fread("chipseq_h3k27ac.300nt.2.bedgraph") %>% setNames(c("chr", "start", "end", "score"))

input.1 <- fread("chipseq_input.300nt.1.1.bedgraph") %>% setNames(c("chr", "start", "end", "score"))
input.2 <- fread("chipseq_input.300nt.2.1.bedgraph") %>% setNames(c("chr", "start", "end", "score"))

h3k27 <- h3k27.1 %>% mutate(score = score + h3k27.2$score)
input <- input.1 %>% mutate(score = score + input.2$score)

chip.h3k27 <- h3k27 %>% mutate(score = log2(score / input$score)) %>% filter(is.finite(score))

write.table(chip.h3k27, "h3k27_chipseq_norm.sum.bedgraph", sep = "\t", col.names = F, row.names = F, quote = F)

