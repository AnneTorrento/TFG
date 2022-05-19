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


print("ITERATING")
plot1 <- list()
plot2 <- list()

for (i in 1:length(obj.list)){
  md <- obj.list[[i]]@meta.data %>% as.data.table
  if(nrow(md) > 150){
    ctype <- unique(md$celltype)
    obj.list[[i]] <- RunPCA(obj.list[[i]], features = features)
    
    obj.list[[i]] <- FindNeighbors(obj.list[[i]], dims = 1:9)
    obj.list[[i]] <- FindClusters(obj.list[[i]], resolution = 0.5)

    #obj.list[[i]] <- RunUMAP(obj.list[[i]], features = features, n.neighbors = 15, verbose = F) 
    obj.list[[i]] <- RunUMAP(obj.list[[i]], features = features, n.neighbors = 70, verbose = F) 
    #INCREASE NUMBER OF NEIGHBORS!!!
    #BUT SOME GROUPS HAVE LESS THAN 17 CELLS, SO MAYBE REMOVE GROUPS WITH LESS THAN 50 CELLS AND SET N.NEIGHBORS TO 50.
    
    plot1[[i]] <-  DimPlot(obj.list[[i]], reduction = "umap") + ggtitle(ctype)
    plot2[[i]] <-  DimPlot(obj.list[[i]], reduction = "umap", group.by = "stage") + ggtitle(ctype)
  }
}


print("CREATING PDF 1")
pdf(file = "UMAP_70nn_Celltypes_group.pdf", 8, 130)
CombinePlots(plot1, ncol = 1)
dev.off()

print("CREATING PDF 2")
pdf(file = "UMAP_70nn_Celltypes_stage.pdf", 8, 130)
CombinePlots(plot2, ncol = 1)
dev.off()


