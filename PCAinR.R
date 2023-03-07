# PCA has 5 steps:
# data normalization, covariance matrix multiplication,
# selection of principal components, Data transformation in new space

#DATA NORMALIZATION

#data normalization ensures that each attribute has the same level
#of contribution, preventing one variable from dominating others.
#For each  variable, normalization is done by subtracting its mean and dividing by its
#standard deviation

#COVARIANCE MATRIX
# this step is computing the covariable matrix from the normzalied data 
#this is a symmetric matrix and each element (i, j) corresponds to the covariance 
#between variables i and j

#EIGENVECTORS AND EIGENVALUES
#an eigenvector represents a direction such as a "vertical" or "90 degrees".
#and eigen value is a number representing the amount of variance present in the 
# data for a given direction. Each eigenvector has its corresponding eigen value


#DATA TRANSFORMATION IN NEW SPACE
#this step re-orients the data onto a new supspace defined by the principal components
#This reorientation is done by multiplying the orginal data by the previously computed 
#eigen vectors 

#R package for correlation analysis. FOcuses on creating and handling R data frames. 
install.packages('corrr')
library('corrr')

#ggcorrplot helps in visualizing the correlation matrix
install.packages('ggcorrplot')
library(ggcorrplot)

#FactoMineR is mainly used for multivariate exploratory data analysis. It gives access to 
#the PCA module to perform prinipal component analysis
install.packages('FactoMineR')
library('FactoMineR')

protein_data <- read.csv("protein.csv")
str(protein_data)

colSums(is.na(protein_data))
numerical_data <- protein_data[,2:10]

head(numerical_data)
data_normalized <- scale(numerical_data)
head(data_normalized)


corr_matrix <- cor(data_normalized)
ggcorrplot(corr_matrix)


data.pca <- princomp(corr_matrix)
summary(data.pca)
data.pca$loadings[, 1:2]


fviz_eig(data.pca, addlabels = TRUE)
install.packages("devtools")
library("devtools")
install_github("kassambara/factoextra")
library("factoextra")

install.packages("scales")
fviz_eig(data.pca, addlabels = TRUE)

fviz_pca_var(data.pca, col.var = "black")
