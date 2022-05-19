.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library("factoextra")
library("FactoMineR")
library(mclust)
library(pvclust)
library(dendextend)
library(Cairo)
library(snow)

# ~~~ upload data ~~~ #
print("READING TABLE")
my_data = read.table(file = "table.txt", header = TRUE, row.names=1) #read the file (raw)
celltypes <- unique(my_data$celltype)
row.names(my_data) <- paste(my_data$sample, my_data$stage, sep="_")

# ~~~ Ward's Hierarchical Clustering ~~~ #
for(i in 1:length(celltypes)){
  print(i)
  print("COMPUTING RESULT")
  result <- pvclust(t(my_data[which(my_data$celltype == celltypes[i]), 1:11]), method.dist="euclidean", method.hclust="ward.D2", nboot=1000, parallel = TRUE)
  
  print("SAVING RDS")
  saveRDS(result, file = paste("pvclust_", celltypes[i],".rds", sep = ""))
}
