library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
head(pbmc)
pbmc.data[1:500,1:30]
pbmc.data[c("CD3D", "RPL11", "MS4A1", "MT-ND1","MT-ND4"), 1:100]
pbmc.data[1:30]
CD3D.data <- pbmc.data[c("CD3D"), 1:2700]
CD3D.data[1:30]
pbmc.data[c("CD3D"), 1:30]
RPL11.data <-- pbmc.data[c("RPL11"), 1:2700]
MS4A1.data <-- pbmc.data[c("MS4A1"), 1:2700]
ND1.data <-- pbmc.data[c("MT-ND1"), 1:2700]
ND4.data <-- pbmc.data[c("MT-ND4"), 1:2700]

data <- data.frame(x = c("CD3D", "RPL11", "MS4A1", 
                         "MT-ND1","MT-ND4"), y= c(CD3D.data, RPL11.data, MS4A1.data, 
                                                  ND1.data, ND4.data))
ggplot(data, aes(x, y)) + 
  geom_point(shape = 21, size = 4, fill = "red") +
  xlab("X-axis Label") +
  ylab("Y-axis Label") +
  ggtitle("Dot Plot")

x_axis_vals <- 1:100
CD3D.data <-- pbmc.data[c("CD3D"), 1:2700]
CD3D.data

library(matrixStats)
rowCounts(as.matrix(CD3D.data))
dim(as.matrix(CD3D.data))
hist_1 <- hist(CD3D.data, xlim=c(6,100), breaks = 40)
hist_2 <- hist(CD3D.data, xlim=c(6, 100), breaks = 40)
hist_2
hist_3 <- hist(CD3D.data, xlim=c(9, 48), ylim = c(-1, 12),breaks = 40)
hist_3
