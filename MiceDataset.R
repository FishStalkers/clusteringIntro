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

mice_obj <- CreateSeuratObject(counts = mm.data, min.cells = 3, min.genes = 200, project = "mm")

head(mice_obj@meta.data)
#everything above this line getting data from google drive
dense.size <- object.size(as.matrix(mm.data))

#mm.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size
head(mice_obj)
mm.data[, 1:30]
head(pbmc.data[1:3, 1:30])
mm.data[c("MIR1302-10", "FAM138A", "OR4F5"), 1:30]
