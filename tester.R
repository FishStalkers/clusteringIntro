library(Seurat)
library(dplyr)
library(Matrix)
library(googledrive)

temp <- tempfile(fileext = "neurons_900_raw_gene_bc_matrices.tar.gz")
dl <- drive_download(
  as_id("https://drive.google.com/file/d/1svEzdxEwh_JuTbsb482odJN6Be6YZSv9/view"), path = temp, overwrite = TRUE)
untar(temp, exdir = tempdir())

#get barcodes
#barcodes <- read.table("C:/Users/Adam/AppData/Local/Temp/RtmpEPZQpX/raw_gene_bc_matrices/mm10/barcodes.tsv")
#head(barcodes)

mm.data <- Read10X(data.dir = "C:/Users/Adam/AppData/Local/Temp/RtmpEPZQpX/raw_gene_bc_matrices/mm10")