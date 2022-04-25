library(Seurat)
library(ggplot2)
library(data.table)
library(magrittr)
library(cytobox)

print("READING SCALED ATLAS")
obj <- readRDS("obj.rds")

print("READING FEATURES")
features <- readRDS("features_LongPolyA.rds")

print("SPLITTING OBJECT BY CELLTYPE")
obj.list <- SplitObject(obj, split.by = "celltype")

print("INITIALIZING VARIABLES")
stages <- unique(obj@meta.data$stage)
celltypes <- unique(obj@meta.data$celltype)

print("ITERATING")
plot1 <- list()

for (i in 1:length(obj.list)){
  md <- obj.list[[i]]@meta.data %>% as.data.table
  if(nrow(md) > 150){    
    obj.list[[i]] <- RunPCA(obj.list[[i]], features = features)

    obj.list[[i]] <- FindNeighbors(obj.list[[i]], dims = 1:2)
    obj.list[[i]] <- FindClusters(obj.list[[i]], resolution = 0.1)
    
    obj.list[[i]] <- RunTSNE(obj.list[[i]], dims = 1:2, perplexity = 50)
    
    plot1[[i]] <-  tsneByMeanMarkerExpression(obj, features, reduction = "tsne", title = celltypes[i], hide_axes = FALSE)
  }
}

print("CREATING PDF")
pdf(file = "tSNE_exprs.pdf", 8, 130)
CombinePlots(plot1, ncol = 1)
dev.off()
