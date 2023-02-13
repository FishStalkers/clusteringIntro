library(Seurat)
library(dplyr)
mm.data <- Read10X(data.dir = "C:/Users/Adam/AppData/Local/Temp/RtmpEPZQpX/raw_gene_bc_matrices/mm10")
mice_obj <- CreateSeuratObject(counts = mm.data, min.cells = 3, min.genes = 200, project = "mm")
mice_obj_vf <- FindVariableFeatures(mice_obj, selection.method ="vst", nfeatures = 200)
head(mice_obj_vf)
top_10_feats <- head(VariableFeatures(mice_obj_vf), 10)
head(mice_obj)


plot1 <- VariableFeaturePlot(mice_obj_vf) + theme(legend.text = element_text(size = 6))
plot2 <- LabelPoints(plot = plot1, points = top_10_feats, repel = TRUE) + theme(legend.text = element_text(size = 6))
#plot1 + plot2
plot2


