#Analysis of Gene Count
pacman::p_load(dplyr, Seurat, patchwork, data.table)
pbmc.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = pbmc.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
head(mice_data)
mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = 
                                                    "^mt-")
dim(mice_data)
expr_raw <- GetAssayData(objeg)
mice_data[rownames(mice_data)]
dataValue <- as.data.frame(as.matrix(mice_data@scale.data)) 
