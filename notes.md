# Introduction 
## Shared Nearest Neighbor
- Traditional clusterring algorithms such as K-means provide a polynomial runtime algorithm to order a set of unlabelled data and find structural similarity (2)
- With higher dimension, the flaws of these traditional algorithms are clear due to the high dimensionality data it may be harder to capture similarities within data in a meaningful manner (2)
  -  The $L_2$ norm is not a good metric to measure similarity of data for high dimensional data
- SNN aims to alleviate this issue:
  - Combines Jarvis-Patrick with DB-scan 
  - Before we delve into the details of SNN, let's first explore the nature of Jarvis-Patrick:
    - Represent the dataset as weighted complete graph $G = (V(G),E(G))$ where edge weight denotes Euclidean distance 
    - Find k-nearest neighbors of all points, and define a function $f(x): V(G) \rightarrow V(G)$ that provides a mapping of some vertex $u$ to $v$ where v denotes the k-nearest neighbor of u 
    - Let's define the setting of mappings as set of tuples $N$
    - We then build our similarity table $S$ based on $N$, let $s_{v,u} \in \R, u\in V(G)$ denote the similarity score of each vertex pair where similarity score is the number of shared neighbors and 0 if k-nearest neighbor isn't shared
    - Once our table S is constructed we can cluster based on the number of shared neighbors and the k parameter
- SNN is a variation of this algorithm that also utilizes DB-scan
- In essence, we precompute the similarity matrix, then we reduce the size of similarity matrix by removing every instance of matrix where similarity is 0
- We then build a graph of the remaining similarity table and select core points to form our neighborhood based on the constraint such that for the KNN list that $sim(u,v) \geq \epsilon $
- We then assign points based on if density is greater than the minimum points and form clusters based on $\epsilon$
- We then align points to cluster
- The overall time complexity of this approach is $O(n^2)$
## Dimensionality Reduction 
### Feature Selection for Individual Datasets: VST
- The following algorithm aims to find differences between features by aiming to find highly variable features
- Simply choosing genes based on the variances between the cells fails to account for the mean variance relationship that is inherent to single-cell RNA-seq (1)
- To compute the mean variance relationship of the data, we can compute the mean and variance of each gene utilizing unnormalized data and applie a $log_{10}$ transformations to both (1)
- Then utilzing regression, we fit a polynomial to predict the variance of each gene or feature as a function of it's mean let's denote this as $f(x): \R \rightarrow \R$ (1)
- Once our regressor is made, we perform normalization based on the data 
- Given a count matrix, let's define $X$ as the row vector of features and $X_{\mu}$ denote the vector corresponding to each mean of the feature
- Let $Z$ define the normalized vector of $X$, moreover let's define a function: $g(X): \R^{|X|} \rightarrow \R^{|X|}, x \in \R^{|X|}$  where our function $g(X)$ performs $f(x)$ on each entry of X. (1)
- We then calculate the standarized variance defined as $\frac{X - X_{\mu}}{g(X)}$ and select features based on this variance (1)
Relevant Sources:
1. https://www.sciencedirect.com/science/article/pii/S0092867419305598?via%3Dihub
2. https://ieeexplore.ieee.org/document/7839671
