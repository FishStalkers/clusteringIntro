
pacman::p_load(dplyr, Seurat, patchwork,bnlearn,Rgraphviz, gRain,sctree)
mice.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = mice.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
val = mice_data[c("Hba-a1","Hba-a2", "Hbb-bt" ,"Hbb-bs" ,"Reln" ,  "Alas2" , 
                  "Ube2c" , "Fabp7" , "Tac2" ,  "Cks2"  )]
df <- as.data.frame(mice_data@assays$RNA@counts
res <- hc(df)
plot(res)

