.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(pvclust)


##############################################################################################################
print("READING DATA")
obj <- readRDS("obj_normalized.rds")
mtGenes <- readRDS("features_LongPolyA.rds")
#nucGenes <- readRDS("nuclearGenes.rds") #list of nuclear genes related to mt expression (known beforehand)
celltypes <- readRDS("celltypes.rds")



##############################################################################################################
print("ADDING HIERARCHICAL (PVCLUST) CLUSTERS TO OBJECT METADATA")
nClusters <- c(5, 6, 4, 4, 5,
               4, 5, 5, 4, 5,
               5, 5, 7, 5, 4,
               5, 7, 4, 7, 5,
               4, 5, 8, 5, 4,
               5, 4, 5, 7, 4,
               5, 4, 4, 2, 4, 
               2, 2)

cluster <- data.frame()

for (i in 1:length(celltypes)){
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- as.data.frame(cutree(result$hclust, k=nClusters[i]))
  rownames(tmp) <- NULL
  cluster <- rbind(cluster, tmp)
}

colnames(cluster) <- "cluster"

obj@meta.data <- cbind(obj@meta.data, cluster)
names(obj@meta.data[ncol(obj@meta.data)]) <- colnames(cluster) #or <- "cluster"

saveRDS(obj, file = "obj_normalized_clusters.rds")



##############################################################################################################
obj <- readRDS("obj_normalized_clusters.rds")

for (i in 1:length(celltypes)){
    #for each stage in each cell type:
        #for each cluster within every stage in a cell type:
            #find differential expression markers
}
