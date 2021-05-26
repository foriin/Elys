library(dm3)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(openxlsx)
s2.chip.bed <- import.bed("/home/artem/IMG/Projects/Elys/chip_capelson/s2_elys_peaks.bed") %>% 
  GRanges()

s2.chip.bg <- import.bedGraph("/home/artem/IMG/Projects/Elys/chip_capelson/s2_elys_chip.bg") %>% 
  GRanges()

kc.nup.npc <- read.xlsx("~/Downloads/Nup98 and HP1a domains in Kc (from Ilyin et al. NAR 2017).xlsx", sheet = 1) %>% 
  select(1:3) %>% setNames(c("chr", "start", "end")) %>% 
  mutate(chr = sub("chr", "", chr), start = start + 1) %>% makeGRangesFromDataFrame()

kc.nup.nuc <- read.xlsx("~/Downloads/Nup98 and HP1a domains in Kc (from Ilyin et al. NAR 2017).xlsx", sheet = 2) %>% 
  select(1:3) %>% setNames(c("chr", "start", "end")) %>% 
  mutate(chr = sub("chr", "", chr), start = start + 1) %>% makeGRangesFromDataFrame()

save(s2.chip.bed, s2.chip.bg, file = "RData/s2.elys.chip.RData")
save(kc.nup.npc, kc.nup.nuc, file = "RData/kc.nup98.domains.RData")

