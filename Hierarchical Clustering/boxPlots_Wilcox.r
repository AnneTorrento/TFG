library(pvclust)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

print("READING TABLE")
my_data = read.table(file = "table.txt", header = TRUE, row.names=1) #read the file (raw)
celltypes <- readRDS("celltypes.rds")

#JITTER BOXPLOTS WITH PAIRWISE COMPARISONS
for(i in 1:length(celltypes)){
  print("LOADING DATA")
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- my_data[which(my_data$celltype == celltypes[i]),]
  tmp$cluster <- cutree(result$hclust, k=7)
  tmp$AvgExprs <- rowMeans(tmp[,2:11])

  print("GENERATING BOXPLOT")
  p <- ggplot(tmp, aes(x=as.factor(cluster), y=AvgExprs)) + 
        geom_boxplot() +
        geom_jitter(shape=16, position=position_jitter(0.2), size = 0.8, aes(col = stage))

  stat.test <- tmp %>% wilcox_test(AvgExprs ~ cluster)
  stat.test <- stat.test %>% add_xy_position(x = "cluster")
  p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)

  pdf(file = paste("WILCOX", celltypes[i], ".pdf", sep = ""))          
  plot(p)
  dev.off()
}
