#Analysis of Gene Count
pacman::p_load(dplyr, Seurat, patchwork, data.table, psych, tidyverse)
mice.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = pbmc.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
#visualizes sum
sumsArray <- rowSums(mice.data)
meansArray <- rowMeans(mice.data)
barplot(sumsArray[1:50], main = "Average of Genes")
barplot(sumsArray[51:100] ,main = "Average of Genes")
barplot(sumsArray[101:150] , main ="Average of Genes")
barplot(sumsArray[151:200], main = "Average of Genes")
barplot(sumsArray[201:250,])
#visualizes means
barplot(meansArray[1:50], main = "Average of Genes")
barplot(meansArray[51:100], main = "Average of Genes")
barplot(meansArray[101:150], main = "Average of Genes")
barplot(meansArray[151:200], main = "Average of Genes")
barplot(meansArray[201:250,], main ="Average of Genes")


head(rownames(mice_data))
head(colnames(mice_data))
#converts to csv and get count of each gene
mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = "^mt-")
mice_datacsv <- GetAssayData(object = mice_data, slot ='counts')
fwrite(as.data.table(sumsArray,keep.rownames = "feature"), "sums.csv")
fwrite(as.data.table(meansArray,keep.rownames = "feature"), "means.csv")
fwrite(as.data.table(mice_datacsv,keep.rownames = "feature"), "counts.csv")
#imports data frame 
df <- read.csv("counts.csv")
describe(df)
dim(df)



