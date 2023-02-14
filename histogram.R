pacman::p_load(dplyr, Seurat, patchwork)
df = read.csv("test.csv")
dim(df)

gene1 <- unlist(df[1,])
hist(as.numeric(unlist(gene1)), main = "Xklr", xlab = "Appearances", ylab = "Count")
gene1 <- df[2,]
hist(as.numeric(unlist(gene1)), main = "Gm1992", xlab = "Appearances", ylab = "Count")
gene1 <- df[3,]
hist(as.numeric(unlist(gene1)), main = "Gm37381", xlab = "Appearances", ylab = "Count")
gene1 <- df[4,]
hist(as.numeric(unlist(gene1)), main = "Rp1", xlab = "Appearances", ylab = "Count")
gene1 <- df[5,]
hist(as.numeric(unlist(gene1)), main = "Rp1.1", xlab = "Appearances", ylab = "Count")
