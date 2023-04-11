
pacman::p_load(dplyr, Seurat, patchwork,bnlearn,Rgraphviz, gRain,sctree)
mice.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = mice.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)

genesVar = c("Hba-a1","Hba-a2", "Hbb-bt" ,"Hbb-bs" ,"Reln" ,  "Alas2" , 
          "Ube2c" , "Fabp7" , "Tac2" ,  "Cks2")
mice_data <- mice_data[genesVar]
print(dim(mice_data))
df <- as.data.frame(mice_data, genes = genesVar)
print(colnames(df))
df <-df[1:10]
res <- hc(df)
plot(res)

