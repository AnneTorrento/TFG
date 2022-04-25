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

print("ITERATING")
plot1 <- list()
plot2 <- list()

for (i in 1:length(obj.list)){
  md <- obj.list[[i]]@meta.data %>% as.data.table
  if(nrow(md) > 150){
    ctype <- unique(md$celltype)
    obj.list[[i]] <- RunPCA(obj.list[[i]], features = features)

    obj.list[[i]] <- FindNeighbors(obj.list[[i]], dims = 1:2)
    obj.list[[i]] <- FindClusters(obj.list[[i]], resolution = 0.1)
   
    obj.list[[i]] <- RunTSNE(obj.list[[i]], dims = 1:2, perplexity = 50)
  
    plot1[[i]] <-  DimPlot(obj.list[[i]], reduction = "tsne") + ggtitle(ctype)
    plot2[[i]] <-  DimPlot(obj.list[[i]], reduction = "tsne", group.by = "stage") + ggtitle(ctype)
  }
}


print("CREATING PDF 1")
pdf(file = "TSNE_group.pdf", 8, 130)
CombinePlots(plot1, ncol = 1)
dev.off()

print("CREATING PDF 2")
pdf(file = "TSNE_stage.pdf", 8, 130)
CombinePlots(plot2, ncol = 1)
dev.off()
