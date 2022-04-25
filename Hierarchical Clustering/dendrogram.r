library("factoextra")
library("FactoMineR")
library(mclust)
library(pvclust)
library(dendextend)
library(Cairo)
library(snow)
library(colourvalues)

# ~~~ upload data ~~~ #
print("READING TABLE")
my_data = read.table(file = "table.txt", header = TRUE, row.names=1) #read the file (raw)
celltypes <- readRDS("celltypes.rds")

# ~~~ Ward's Hierarchical Clustering ~~~ #
for(i in 1:length(celltypes)){
  print("LOADING RESULT")
  result <- readRDS(paste("HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))

  print("FORMING DENDROGRAM")
  dend = as.dendrogram(result)
  
  colors_to_use <- as.numeric(my_data[which(my_data$celltype == celltypes[i]),]$stage)
  colors_to_use <- colors_to_use[order.dendrogram(dend)]
  colors_to_use <- colour_values(colors_to_use)
  labels_colors(dend) <- colors_to_use

  print("CREATING SVG")
  svg(filename = paste("dend_", celltypes[i], ".svg", sep = ""), width=50)
  plot(dend, hang=-1, cex=0.5)
  dev.off()
}
