setwd("C:/Users/hrmin/Desktop/VIP2023")
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

download.file("https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", 'pbmc3k')
untar('pbmc3k', files = 'filtered_gene_bc_matrices/')
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
head(pbmc)
pbmc.data[1:500,1:30]
pbmc.data[c("CD3D", "RPL11", "MS4A1", "MT-ND1","MT-ND4"), 1:100]

CD3D.data <- pbmc.data[c("CD3D"), 1:2700]
CD3D.data
RPL11.data <- pbmc.data[c("RPL11"), 1:2700]
MS4A1.data <- pbmc.data[c("MS4A1"), 1:2700]
ND1.data <- pbmc.data[c("MT-ND1"), 1:2700]
ND4.data <- pbmc.data[c("MT-ND4"), 1:2700]

data <- data.frame(x = c("CD3D", "RPL11", "MS4A1", 
                         "MT-ND1","MT-ND4"), y= c(CD3D.data, RPL11.data, MS4A1.data, 
                                                  ND1.data, ND4.data))
ggplot(data, aes(x, y)) + 
  geom_point(shape = 21, size = 4, fill = "red") +
  xlab("X-axis Label") +
  ylab("Y-axis Label") +
  ggtitle("Dot Plot")
hist_1 <- hist(CD3D.data, xlim=c(9, 48), ylim = c(0, 12),breaks = 40)
hist_2 <- hist(RPL11.data, xlim=c(9, 70), ylim = c(0, 225),breaks = 40)
hist_3 <- hist(MS4A1.data, xlim=c(0, 36), ylim = c(0, 100),breaks = 40)
hist_4 <- hist(ND1.data, xlim=c(9, 50), ylim = c(0, 100),breaks = 40)
hist_5 <- hist(ND4.data, xlim=c(10, 40), ylim = c(0, 100),breaks = 30)
hist_1
hist_2
hist_3
hist_4
hist_5