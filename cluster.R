pacman::p_load(dplyr, Seurat, patchwork)
#loads data 
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
head(pbmc.data)
#creates object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 
          3)
#plots feature of cells
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = 
                          "nFeature_RNA")
plot1 + plot2        
