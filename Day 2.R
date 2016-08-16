library(curl)
my_data <- read.table(curl("https://tcga-data.nci.nih.gov/docs/publications/gbm_exp/unifiedScaledFiltered.txt"),sep="\t")

##A. Calculate the pairwise correlations between genes. What are the dimensions of the output matrix? What would it be if we correlated samples instead of genes?

w <- cor(t(my_data))
head(w)

##B. What is the range of values in this matrix? What are the highest and lowest values you would expect and why?
dim(w)
max(w) #genes correlated with themselves will have value of 1
min(w)

##C. Find the pair of genes that are most correlated. Do you need to operate over the entire matrix to do this?
w2 <- as.data.frame(as.table(w))
head(w2)
w3 <- w2[w2$Var1 != w2$Var2, ]
head(w3)
w4 <- w3[order(w3$Freq, decreasing = TRUE), ]
head(w4)

##D. By default, the R function cor() uses the Pearson method for calculating correlation coefficients. If we change to using a rank-based method like Spearman’s rho, do we find a different pair of genes that are the most highly correlated?
my_data_t <- t(my_data)
spear_w <- cor(my_data_t, method = "spearman")
spear_w2 <- as.data.frame(as.table(spear_w))
head(spear_w2)
spear_w3 <- spear_w2[spear_w2$Var1 != spear_w2$Var2, ]
head(spear_w3)
spear_w4 <- spear_w3[order(spear_w3$Freq, decreasing = TRUE), ]
head(spear_w4)

##E. Plot a matrix of the top 20 most correlated genes (hint: consider using the library corrplot)
install.packages("corrplot")
library(corrplot)
top20_var1 <- unique(w4$Var1[1:20])
top20_var2 <- unique(w4$Var2[1:20])
result <- w[top20_var1, top20_var2]
corrplot(result)

##Clustering
#A. Cluster the top 20 most correlated genes using complete hierarchical clustering and visualize the output tree
d <- dist(result)
cluster_d<- hclust(d)
plot(cluster_d)

#B. Try using different clustering methods and observe the effects upon how correlated genes are grouped together.

#C. If we cluster these genes with complete hierarchical clustering and divide the tree into 4 groups, how many genes are in the largest group? The smallest?
rect.hclust(cluster_d, k=4, border="red")
