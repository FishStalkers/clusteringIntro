download.file("https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
              , 'pbmc3k')
untar('pbmc3k', files = 'filtered_gene_bc_matrices/')
barcodes <- read.table("filtered_gene_bc_matrices/hg19/barcodes.tsv")
head(barcodes)

genes <- read.delim("filtered_gene_bc_matrices/hg19/genes.tsv", header = FALSE)
head(genes)
library(Matrix)
mat <- readMM(file = "filtered_gene_bc_matrices/hg19/matrix.mtx")
head(mat)
mat[1:5, 1:10]

library(dplyr)
library(Seurat)

#load the PBMC dataset

#Read10X() : This function is from the Seurat package and will use 
#the Cell Ranger output directory as input. In this way individual 
#files do not need to be loaded in, instead the function will load 
#and combine them into a sparse matrix for you.

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19")

#initialize Suerat Object with the raw (non-normalized data)
#The Seurat object serves as a container that contains both data 
#(like the count matrix) and analysis (like PCA, or clustering results) 
#for a single-cell dataset. Before using Seurat to analyze scRNA-seq data, 
#we can first have some basic understanding about the Seurat object from here.

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 5)
#Draws a violin plot of single cell data (gene expression, metrics, PC scores, etc.)
#object, features cols = NULL, ...
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#We can take a look and see that the unique counts and feature are correlated. 
#In addition, we see that low counts appear to correlate with high mitochondrial mapping percentage.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1 + plot2

RidgePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)

# tells us what the "VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)"
# does but on an easier level

#To identify thresholds
# use MAD

library(scater)

qc.nCount_RNA <- isOutlier(pbmc$nCount_RNA, log = TRUE, type = "both")
qc.nFeature_RNA <- isOutlier(pbmc$nFeature_RNA, log = TRUE, type = "both")
qc.percentmt <- isOutlier(pbmc$percent.mt, type = "both")

attr(qc.nCount_RNA, "thresholds")
#head(qc.nCount_RNA)
attr(qc.nFeature_RNA, "thresholds")
attr(qc.percentmt, "thresholds")
#based off the previous attriubutes from min and max. Take out outliers that are more than certain thresholds or less
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
#normalize the data
head(pbmc@meta.data)
summary(pbmc@meta.data)
pbmc <- NormalizeData(pbmc)

#Find Most Variable Features
#It is common to identify highly variable features or genes for 
#dimensional reduction. By reducing your analysis to the highly 
#variable genes, you account for most of the biological heterogeneity 
#or factors in your data and hopefully ignore a majority of the noise 
#while reducing computational work and time. As such, the highly variable 
#genes should enable us to isolate the real biological signals.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
top10
plot1 <- VariableFeaturePlot(pbmc) + theme(legend.text = element_text(size = 6))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.text = element_text(size = 6))
plot2

#scaling the data
#All gene expression will be centered to 0.
#Scales the expression of each gene to have a variance of 1 
#so all genes have equal contributions
nrow(pbmc@meta.data)                                                                        
pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc)) #VariableFeatures is used to call the highly variable genes from the object

#Perform PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
print(pbmc[["pca"]], dims = 1:10, nfeatures = 10)

pltt1 <- DimPlot(pbmc, reduction = "pca")
pltt2 <- ElbowPlot(pbmc)
pltt1 + pltt2


#clustering the cells
#To aid in summarizing the data for easier interpretation, scRNA-seq 
#is often clustered to empirically define groups of cells within the data
#that have similar expression profiles.

#Seurat uses a graph-based clustering approach. There are additional 
#approaches such as k-means clustering or hierarchical clustering.
install.packages('installr')
library(installr)
install.Rtools()
install.packages('Matrix')

library(Matrix)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(pbmc)
#embedd clusters onto 2D Space
pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)
#visualize umap, tsne, and pca
library(ggplot2)
plot1 <- DimPlot(pbmc, reduction = 'tsne')
plot2 <- DimPlot(pbmc, reduction = 'umap')
plot3 <- DimPlot(pbmc, reduction = 'pca')
plot1 +  plot2 + plot3

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

