library(pvclust)
library(ggplot2)
library(gridExtra)

print("READING TABLE")
my_data = read.table(file = "table.txt", header = TRUE, row.names=1)
celltypes <- readRDS("celltypes.rds")

#JITTER BOXPLOTS
for(i in 1:length(celltypes)){
  print("LOADING DATA")
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- my_data[which(my_data$celltype == celltypes[i]),]
  tmp$cluster <- cutree(result$hclust, k=7)
  tmp$AvgExprs <- rowMeans(tmp[,2:11])

  print("GENERATING BOXPLOT")
  jitter <- ggplot(tmp, aes(x=as.factor(cluster), y=AvgExprs)) + 
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 0.8, aes(col = stage))

  pdf(file = paste("HClustBPlot_", celltypes[i], ".pdf", sep = ""))          
  plot(jitter)
  dev.off()
}
