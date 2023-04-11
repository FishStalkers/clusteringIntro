install.packages('SeuratObject')
install.packages("pacman")
library(pacman)
p_load('SeuratObject')  
library(Seurat)
library(bnlearn)
library(parallel)
library(data.table)
library(readxl)
counts <- as.data.frame(as.matrix(bb_demux_022322$RNA@scale.data))

counts <- as.data.frame(counts)


counts.dedup <- dedup(as.data.frame(counts), .95, debug = FALSE)

set.seed(42)
sink(file = "learning-log-hc.boot.strength.txt")

dag.hybrid.group <- hc(counts)
cl <- parallel::makePSOCKcluster(16,outfile="debug.txt")
boot.hc.hybrid.group = boot.strength(data=counts,R=1000,algorithm="hc", algorithm.args = list(cluster=cl))

sink()

stopCluster(cl)


avg.boot.hc.hybrid.group <- averaged.network(boot.hc.hybrid.group, threshold=0.5)
avg.boot.hc.dag.hybrid.group <- cextend(avg.boot.hc.hybrid.group)
fitted.hybrid.group <- bn.fit(avg.boot.hc.dag.hybrid.group,counts)


detach("package:bnlearn")
save.image(file="Network_1000rep_BB.RData")
quit()