library(Seurat)
library(ggplot2)

print("Loading LONG POLY-A TAIL FEATURES")
features <- readRDS("features_LongPolyA.rds")

print("READING SCALED ATLAS")
atlas <- readRDS("ATLAS6_scaled.rds")

print("CHANGING IDENT TO STAGE")
atlas <- SetIdent(atlas, value = atlas@meta.data$stage) 

print("SPLITTING OBJECT BY CELLTYPE")
obj.list <- SplitObject(atlas, split.by = "celltype")

print("ITERATING")
plot1 <- list()
for (i in 1:length(obj.list)){
  plot1[[i]] <-  DoHeatmap(obj.list[[i]], features = features, size = 4, angle = 90) + 
                          NoLegend() + 
                          ggtitle(unique(obj.list[[i]]@meta.data$celltype))
}

print("CREATING PDF")
pdf(file = "Heatmaps.pdf", 8, 130)
CombinePlots(plot1, ncol = 1, legend = "none")
dev.off()
