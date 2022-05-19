.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(ggplot2)

obj <- readRDS("obj_scaled.rds")

print("READING FEATURES")
features <- readRDS("features_LongPolyA.rds")

print("SPLITTING OBJECT BY CELLTYPE")
obj.list <- SplitObject(obj, split.by = "celltype")


print("ITERATING")
plot1 <- list()
plot2 <- list()

for (i in 1:length(obj.list)){
  obj.list[[i]] <- RunPCA(obj.list[[i]], features = features)
  
  pdf(file = paste("Elbow_", i, ".pdf", sep = ""), 8, 7)
  ElbowPlot(obj.list[[i]])
  dev.off()

  #obj.list[[i]]@meta.data <- obj.list[[i]]@meta.data[, -which(colnames(obj.list[[i]]@meta.data) %in% "cluster")]
  #obj.list[[i]]@meta.data <- obj.list[[i]]@meta.data[, -which(colnames(obj.list[[i]]@meta.data) %in% "umapX")]
  #obj.list[[i]]@meta.data <- obj.list[[i]]@meta.data[, -which(colnames(obj.list[[i]]@meta.data) %in% "umapY")]
  
  #print(colnames(obj.list[[i]]@meta.data))
  obj.list[[i]] <- RunUMAP(obj.list[[i]], features = features, n.neighbors = 15, verbose = F) 
  print(colnames(obj.list[[i]]@meta.data))
  
  plot1[[i]] <-  DimPlot(obj.list[[i]], reduction = "umap")
  plot2[[i]] <-  DimPlot(obj.list[[i]], reduction = "umap", group.by = "stage")
}

print("CREATING PDF 1")
pdf(file = "UMAP_Celltypes_group.pdf", 8, 130)
CombinePlots(plot1, ncol = 1)
dev.off()

print("CREATING PDF 2")
pdf(file = "UMAP_Celltypes_stage.pdf", 8, 130)
CombinePlots(plot2, ncol = 1)
dev.off()


