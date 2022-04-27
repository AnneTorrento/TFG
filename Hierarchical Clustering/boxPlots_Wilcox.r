library(pvclust)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(rstatix)

my_data = read.table(file = "table.txt", header = TRUE, row.names=1) #read the file (raw)
celltypes <- readRDS("celltypes.rds")


#############################################################
#JITTER BOXPLOTS WITH PAIRWISE COMPARISONS
#############################################################
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

#############################################################
#JITTER BOXPLOTS IN GRID PER mt-GENE
#############################################################
library(devtools)
library(gridExtra)

features <- readRDS("features_LongPolyA.rds")
features <- gsub("-", ".", features)

for(i in 1:length(celltypes)){
  print("LOADING DATA")
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- my_data[which(my_data$celltype == celltypes[i]),]
  tmp$cluster <- cutree(result$hclust, k=7)
  plots <- list()

  for(j in 1:length(features)){
    gene <- as.character(features[j])

    print("GENERATING BOXPLOT")

    p1 <- eval(substitute(
      ggplot(tmp, aes(x = as.factor(cluster), y = get(colnames(tmp)[1+j]))) +
      geom_boxplot() +
      geom_jitter(shape=16, position=position_jitter(0.2), size = 0.8, aes(col = stage)) +
      labs(title = gene) +
      xlab("Cluster") +
      ylab("Expression")
      , list(j = j)))

    plots[[j]] <- p1
  }

  print("COMPILING PDF")
  p <- grid.arrange(grobs = plots, ncol = 2)
  pdf(file = paste("expression_", celltypes[i], ".pdf", sep = ""), 8, 20)
  plot(p)
  dev.off()
}



#############################################################
#BOXPLOTS WITH PAIRWISE COMPARISONS IN A GRID BY STAGE
#############################################################
for(i in 1:length(celltypes)){
  print("LOADING DATA")
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- my_data[which(my_data$celltype == celltypes[i]),]
  tmp$cluster <- cutree(result$hclust, k=7)
  tmp$AvgExprs <- rowMeans(tmp[,2:11])
  plots <- list()
  
  for(j in 1:length(unique(tmp$stage))){
    stage <- which(tmp$stage == unique(tmp$stage)[j])
    
    if(stage){
      tmp_stage <- as.data.frame(cbind(tmp[stage,15],
                                       tmp[stage,13],
                                       tmp[stage,14]))
      colnames(tmp_stage) <- c("AvgExprs", "stage", "cluster")
      tmp_stage <- transform(tmp_stage, AvgExprs = as.numeric(AvgExprs))
      
      print("GENERATING BOXPLOT")
      p1 <- eval(substitute(
        ggplot(tmp_stage, aes(x = as.factor(cluster), y = AvgExprs)) +
          geom_boxplot(aes(fill = cluster)) +
          labs(title = tmp_stage$stage[1]) +
          xlab("Cluster") +
          ylim(0, 4)
        , list(j = j)))
      
      if(length(unique(tmp_stage$cluster)) == 2){
        comparisons <- list(unique(tmp_stage$cluster))
        p1 <- p1 + geom_signif(comparisons = comparisons, 
                               map_signif_level=TRUE)
      }
      else if(length(unique(tmp_stage$cluster)) > 2){
        stat.test <- tmp_stage %>% wilcox_test(AvgExprs ~ cluster)
        stat.test <- stat.test %>% add_xy_position(x = "cluster")
        p1 <- p1 + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)  
      }
      
      plots[[j]] <- p1
    }
  }
  
  print("COMPILING PDF")
  cols <- length(plots)
  p <- grid.arrange(grobs = plots, ncol = cols)
  pdf(file = paste("WILCOX_grid_", celltypes[i], ".pdf", sep = ""), 20, 4)
  plot(p)
  dev.off()
}
