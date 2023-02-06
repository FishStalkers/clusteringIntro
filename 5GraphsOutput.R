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

CD3D.data <-- pbmc.data[c("CD3D"), 1:2700]
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
