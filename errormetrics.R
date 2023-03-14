library(Seurat)
library(dplyr)
library(Matrix)
library(googledrive)
library(dplyr)


temp <- tempfile(fileext = "neurons_900_raw_gene_bc_matrices.tar.gz")
dl <- drive_download(
  as_id("https://drive.google.com/file/d/1svEzdxEwh_JuTbsb482odJN6Be6YZSv9/view"), path = temp, overwrite = TRUE)
untar(temp, exdir = tempdir())

mm.data <- Read10X(data.dir = "C:/Users/hrmin/Desktop/VIP2023/raw_gene_bc_matrices/mm10")

head(mm.data[2,1:5])
mice_obj <- CreateSeuratObject(counts = mm.data, min.cells = 3, min.genes = 200, project = "mm")

head(mice_obj@meta.data)

#potential error metrics to use

observed <- c(1, 3, 5, 7, 9)
predicted <- c(2, 3, 6, 8, 10)

df <- data.frame(observed, predicted)
# mean absolute error 
mae <- df %>%
  summarise(MAE = mean(abs(observed - predicted)))
  
# mean squared error
mse <- df %>%
  summarise(MSE = mean((observed - predicted)^2))

#mean squared error
mse <- df %>%
  summarise(MSE = mean((observed - predicted)^2))