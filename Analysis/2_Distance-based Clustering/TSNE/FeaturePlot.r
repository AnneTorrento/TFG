.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(ggplot2)
library(data.table)
library(magrittr)

print("READING SCALED ATLAS")
obj <- readRDS("obj.rds")

print("READING FEATURES")
features <- readRDS("features_LongPolyA.rds")

print("SPLITTING OBJECT BY CELLTYPE")
obj.list <- SplitObject(obj, split.by = "celltype")

# print("SPLITTING OBJECT BY STAGE")
# obj.list <- SplitObject(obj, split.by = "stage")

print("INITIALIZING VARIABLES")
stages <- unique(obj@meta.data$stage)
celltypes <- unique(obj@meta.data$celltype)


print("ITERATING")
for (i in 1:length(celltypes)){
  for (j in 1:length(stages)){    
    plots <- list()
    obj.list[[i]] <- SetIdent(obj.list[[i]], value = obj.list[[i]]@meta.data$stage)
    if (celltypes[i] %in% obj.list[[i]]@meta.data$celltype){
      ob <- subset(obj.list[[i]], idents = celltypes[i])
    
      md <- ob@meta.data %>% as.data.table
      if(nrow(md) > 150){
        ob <- RunPCA(ob, features = features)

        ob <- FindNeighbors(ob, dims = 1:2)
        ob <- FindClusters(ob, resolution = 0.1)
   
        ob <- RunTSNE(ob, dims = 1:2, perplexity = 50)

        plots[[j]] <- FeaturePlot(ob, features = features) + ggtitle(celltypes[i])
      }
    }
  }
  if(celltypes[i] == "Forebrain/Midbrain/Hindbrain"){
    pdf(file = "FP_Forebrain_Midbrain_Hindbrain.pdf", 20, 130)  
    CombinePlots(plots, ncol = 1)
    dev.off()
  }
  else{
    pdf(file = paste("FP_", celltypes[i], ".pdf", sep = ""), 20, 130)  
    CombinePlots(plots, ncol = 1)      
    dev.off()
  }
}


# print("ITERATING")
# for (i in 1:length(stages)){
#   for (j in 1:length(celltypes)){    
#     obj.list[[i]] <- SetIdent(obj.list[[i]], value = obj.list[[i]]@meta.data$celltype)
#     if (celltypes[j] %in% obj.list[[i]]@meta.data$celltype){
#       ob <- subset(obj.list[[i]], idents = celltypes[j])
    
#       md <- ob@meta.data %>% as.data.table
#       if(nrow(md) > 150){
#         ob <- RunPCA(ob, features = features)

#         ob <- FindNeighbors(ob, dims = 1:2)
#         ob <- FindClusters(ob, resolution = 0.1)
   
#         ob <- RunTSNE(ob, dims = 1:2, perplexity = 50)
#         if(celltypes[j] == "Forebrain/Midbrain/Hindbrain"){
#           pdf(file = paste("FP_TSNE_Forebrain_Midbrain_Hindbrain_", stages[i], ".pdf", sep = ""), 20, 10)
#           print(FeaturePlot(ob, features = features))
#           dev.off()
#         }
#         else{
#           pdf(file = paste("FP_TSNE_", celltypes[j], "_", stages[i], ".pdf", sep = ""), 20, 10)
#           print(FeaturePlot(ob, features = features))
#           dev.off()
#         }
#       }
#     }
#   }
# }
