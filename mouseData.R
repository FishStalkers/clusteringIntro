pacman::p_load(dplyr, Seurat, patchwork)
#reads data 
pbmc.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = pbmc.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
head(mice_data)
        mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = 
                                                            "^mt-")
#plotting out features for filtering
VlnPlot(mice_data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
        ncol = 3)
#plot out distribution of data
plot1 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = 
                          "percent.mt", pts.size = 0.5)
plot2 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = 
                          "nFeature_RNA",pts.size = 0.5,  plot.cor = TRUE)
plot1 + plot2
#filter out outliars
mice_data <- subset(mice_data, subset = nFeature_RNA > 200 & nCount_RNA < 25000 & percent.mt < 5)
mice_data <- NormalizeData(mice_data, normalization.method = "LogNormalize")
mice_data <- FindVariableFeatures(mice_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mice_data),10)
top10
#finds variable features
plot1 <- VariableFeaturePlot(mice_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scales and normalize data
all.genes <- rownames(mice_data)
mice_data <- ScaleData(mice_data)
#reduce data 
mice_data <- RunPCA(mice_data, features = VariableFeatures(object = mice_data))
print(mice_data[["pca"]], dims = 1:2, nfeatures = 14)

#visualizes impact of each gene
VizDimLoadings(mice_data, dims = 1:2, nfeatures = 14, reduction = "pca")
#dim plot
DimPlot(mice_data, reduction ="pca")
DimHeatmap(mice_data, dims = 1, cells = 500, balanced = TRUE)
#elbow method 
ElbowPlot(mice_data)
#Clusterring
mice_data <- FindNeighbors(mice_data, dims = 1:13)
mice_data <- FindClusters(mice_data, resolution = 0.5)
#UMAP
mice_data <- RunUMAP(mice_data, dims = 1:10)
DimPlot(mice_data, reduction = "umap")
#Finding Markers:
mice_data.markers <- FindAllMarkers(mice_data, only.pos =  TRUE, min.pct = 0.25, logfc.threshold = 0.25)
x <- mice_data.markers %>% group_by(cluster) %>% top_n(n=1, wt = avg_log2FC)
plot1 <- FeaturePlot(mice_data, features = x$gene[1:4])
plot2 <- FeaturePlot(mice_data, features = x$gene[5:8])
plot1 + plot2
top10 <- mice_data.markers %>% group_by(cluster) %>% top_n(n=1, wt = avg_log2FC)
top10
