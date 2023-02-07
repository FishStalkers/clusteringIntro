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
plot1 <- VariableFeaturePlot(mice_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(mice_data)
