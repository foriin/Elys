library(reshape2)
library(stringr)
library(dm3)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(openxlsx)


s.genes.npc <- read.xlsx("tables/Genes overlapped with only Elys_NPC or Elys_nucl and with their combinations.xlsx",
                     sheet = 2) %>% filter(!is.na(start)) %>%  makeGRangesFromDataFrame()
a.genes.npc <- read.xlsx("tables/Genes overlapped with only Elys_NPC or Elys_nucl and with their combinations.xlsx", 
                     sheet = 3) %>% filter(!is.na(start)) %>% makeGRangesFromDataFrame()
export.bed(s.genes.npc, "bed/silent.genes.elys.x.npc.bed")
export.bed(a.genes.npc, "bed/active.genes.elys.x.npc.bed")

dt.tabs <- lapply(dir("metagene/silact/", pattern = "tab$"), function(file){
  tab <- read.table(file.path("metagene/silact", file), skip = 1)[, -1] %>% t()
  
})
names(dt.tabs) <- gsub("metaplot\\.(.*)elys.npc\\.(.+)\\.tab", "\\1\\2",
                       dir("metagene/silact", pattern = "tab$"))
names(dt.tabs)[grepl("(active|silent).gaf.300", names(dt.tabs))] <- sub("(active|silent)(.gaf.300)", "\\1.emb\\2", names(dt.tabs)[grepl("(active|silent).gaf.300", names(dt.tabs))])


dt.tabs.2 <- lapply(unique(sub("^active\\.|silent\\.(.+)$", "\\1", names(dt.tabs))),
                    function(idx){
                      pr <- as.data.frame(do.call(cbind,
                                                  dt.tabs[grepl(idx,
                                                                names(dt.tabs))])) %>%
                        setNames(c("active", "silent")) %>%
                        mutate(kb = 1:600)

                      melt(pr, id.vars = "kb")
                    })
names(dt.tabs.2) <- unique(sub("^active\\.|silent\\.(.+)$", "\\1", names(dt.tabs)))

pdf("plots/metaplots.for.elys.npc.silent.or.active.genes.pdf")
for (i in unique(sub("^active\\.|silent\\.(.+)$", "\\1", names(dt.tabs)))){
  p <- ggplot(dt.tabs.2[[i]], aes(x = kb, y = value))+
    geom_line(aes(col = variable))+
    # scale_color_manual(values = c("blue",
    #                               "green",
    #                               "purple"))+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          legend.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_x_continuous(breaks=c(1,50,550,600),
                       labels=c("-500", "start", "end", "+500"))+
    ylab(i)+
    ggtitle("Genes, overlapping with ELYSxNPC, silent or active")
    print(p)
}
dev.off()






