library(rtracklayer)

setwd("~/IMG/Projects/Elys/kd_rnaseq/bed/")

sub("\\.200nt.*$", "\\.nreads\\.mapped\\.unique\\.txt", dir(pattern = "200nt")[1])


for (i in dir(pattern = "200nt")){
  bg <- import.bedGraph(i)
  nr <- scan(sub("\\.200nt.*$", "\\.nreads\\.mapped\\.unique\\.txt", dir(pattern = "200nt")[1]),
             double())
  bg$score <- bg$score / nr * 10e6
  export.bedGraph(bg, sub("bedgraph", "rpm\\.bedgraph", i))
}



ec <- import.bedGraph("emb.cnt.200nt.bedgraph")
ec.nr <- scan("emb.cnt.nreads.mapped.unique.txt", double())
ec$score <- ec$score / ec.nr * 10e6
export.bedGraph(ec, "emb.cnt.200nt.rpm.bedgraph")

ek <- import.bedGraph("Emb_kdelys.200nt.bedgraph")
ek.nr <- scan("emb.kdelys.nreads.mapped.unique.txt", double())
ek$score <- ek$score / ek.nr * 10e6
export.bedGraph(ek, "emb.kdelys.200nt.rpm.bedgraph")

length(ek$score)
length(ec$score)
cor(ek$score, ec$score)


s2c <- import.bedGraph("S2_control.200nt.bed")
s2c.nr <- scan("s2.control.nreads.mapped.unique.txt", double())
s2c$score <- s2c$score / s2c.nr * 10e6
export.bedGraph(s2c, "s2.control.200nt.rpm.bedgraph")

s2k <- import.bedGraph("S2_kdelys.200nt.bedgraph")
s2k.nr <- scan("s2.kdelys.nreads.mapped.unique.txt", double())
s2k$score <- s2k$score / s2k.nr * 10e6
export.bedGraph(s2k, "s2.kdelys.200nt.rpm.bedgraph")

cor(s2k$score, s2c$score)
cor(ek$score, s2k$score)
