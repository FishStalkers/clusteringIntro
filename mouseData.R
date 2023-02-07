pacman::p_load(dplyr, Seurat, patchwork)
#reads data 
pbmc.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = pbmc.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
head(mice_data)
        mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = "^mt-")
#plotting out features for filting
VlnPlot(mice_data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
        ncol = 3)
#plot out distribution of data
plot1 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#filter out outliers
mice_data <- subset(mice_data, subset = nFeature_RNA > 200 & ncount_RNA < 25000 & percent.mt < 5)
