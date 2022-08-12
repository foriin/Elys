library(vioplotx)
library(openxlsx)
library(dplyr)

vol <- read.xlsx("tables/Volume of nuclei.xlsx", startRow = 2)
dist <- read.xlsx("tables/Distances from FISH signals.xlsx", startRow = 2)
dist.nl <- dist[, 1:2]
dist.fish <- dist[, 3:4]

pdf("plots/vioplots.elyskd.pdf")
vioplotx(vol, names = names(vol),
         col = c("green", "magenta"), main = "Nuclei volumes", ylab = expression(paste("Volume, ", mu, "m")^3))
vioplotx(dist.nl, names = names(dist.nl), col = c("green", "magenta"),
         main = "Distances between FISH signals and NL",
         ylab = expression(paste("Distance, ", mu, "m")))
vioplotx(dist.fish, names = names(dist.fish), col = c("green", "magenta"),
         main = "Distances between FISH signals",
         ylab = expression(paste("Distance, ", mu, "m")))
dev.off()

ge.tpm <- read.xlsx("tables/expression for box plots.xlsx", sheet = 1, startRow = 2)


ge.fc <- tibble(read.xlsx("tables/expression for box plots.xlsx", sheet = 2, startRow = 2))


ge.fc.int <- tibble(read.xlsx("tables/expression for box plots.xlsx", sheet = 3, startRow = 2))

pdf("plots/elys.nup98.ge.pdf", height = 8, width = 8)
  boxplot(ge.tpm %>% setNames(sub("\\.", " ", names(ge.tpm))), outline = F,
          ylab = "TPM")
  boxplot(ge.fc %>% setNames(sub("\\.", " ", names(ge.fc))), outline = F,
          ylab = "Fold Change")
  boxplot(ge.fc.int %>% setNames(sub("\\.", " ", names(ge.fc.int))), outline = F,
          ylab = "Fold Change")
dev.off()

