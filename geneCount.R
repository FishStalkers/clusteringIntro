#Analysis of Gene Count
pacman::p_load(dplyr, Seurat, patchwork, data.table)
pbmc.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = pbmc.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
head(rownames(mice_data))
head(colnames(mice_data))
mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = "^mt-")
mice_datacsv <- GetAssayData(object = mice_data, slot ='counts')
fwrite(as.data.table(mice_datacsv, ))
dim(mice_data)
