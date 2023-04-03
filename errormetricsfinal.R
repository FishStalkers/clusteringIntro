#install.packages("caret")
library(MASS)
library(caret)
library(Rcpp)
# Load the Boston Housing dataset
data(Boston)

# Split the dataset into training and test sets
set.seed(123)
trainIndex <- createDataPartition(Boston$medv, p = 0.8, list = FALSE)
training <- Boston[trainIndex,]
testing <- Boston[-trainIndex,]

# Train a linear regression model
model <- train(medv ~ ., data = training, method = "lm")

# Make predictions on the test set
predictions <- predict(model, newdata = testing)

# Compute the principal components of the test set
pca <- prcomp(testing[, -14], scale. = TRUE)

# Project the test set onto the first two principal components
pc <- predict(pca, testing[, -14])[, 1:2]

# Compute the MSE of the linear regression model on the first two principal components
mse <- mean((predict(model, as.data.frame(pc)) - testing$medv)^2)
mse

#-----------------------------------------------------------------
data <- read.csv("data.csv", header = TRUE)
kmeans_model <- kmeans(data, centers = 3)

silhouette(kmeans_model$cluster, dist(data))

#other error metrics for unlabeled data: Dunn index,
#Calinski-Harabasz index, and Davies-Bouldin index
pacman::p_load(dplyr, Seurat, patchwork)
#reads data 
mice.data <- Read10X(data.dir = "raw_gene_bc_matrices/mm10/")
mice_data <- CreateSeuratObject(counts = mice.data, project  ="mice_data", 
                                min.cells = 3, min.features = 200)
head(mice_data)
mice_data[["percent.mt"]] <- PercentageFeatureSet(mice_data, pattern = 
                                                    "^mt-")
head(rownames(mice_data)) 
head(mice_data["Xkr4" ]["nCount_RNA"])
#plotting out features for filtering
VlnPlot(mice_data, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
        ncol = 3)
#plot out distribution of data
plot1 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = 
                          "percent.mt", pts.size = 0.5)
plot2 <- FeatureScatter(mice_data, feature1 = "nCount_RNA", feature2 = 
                          "nFeature_RNA",pts.size = 0.5,  plot.cor = TRUE)
plot1 + plot2
#filter out outliars
mice_data <- subset(mice_data, subset = nFeature_RNA > 200 & nCount_RNA < 25000 & percent.mt < 5)
mice_data <- NormalizeData(mice_data, normalization.method = "LogNormalize")
mice_data <- FindVariableFeatures(mice_data, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(mice_data),10)
top10
#finds variable features
plot1 <- VariableFeaturePlot(mice_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scales and normalize data
all.genes <- rownames(mice_data)
mice_data <- ScaleData(mice_data)
#reduce data 
mice_data <- RunPCA(mice_data, features = VariableFeatures(object = mice_data))
print(mice_data[["pca"]], dims = 1:2, nfeatures = 14)