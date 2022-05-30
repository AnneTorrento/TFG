.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(ggplot2)
library(xlsx)
library(data.table)
library(gridExtra)
suppressPackageStartupMessages({
  library(rlang)
})

##############################################################################################################
print("READING DATA")

obj <- readRDS("obj_scaled_clusters.rds")
signifClust <- readRDS("signifDiffClusters.rds")
all_celltypes <- readRDS("celltypes.rds")
celltypes <- unique(signifClust$celltype)
mtGenes <- readRDS("features_LongPolyA.rds")



##############################################################################################################
print("CREATING HEATMAPS")

DE_counter <- 0
DE_features <- data.frame(celltype = character(), clusters = character(), feature = character())

index <- 1
obj.list <- SplitObject(obj, split.by = "celltype") # SPLITTING OBJECT BY CELLTYPE

for (i in 1:length(all_celltypes)){
  if (all_celltypes[i] != celltypes[index]){
    next
  }

  celltype <- celltypes[index]
  obtmp <- obj.list[[i]]
  tmpdf <- signifClust[which(signifClust$celltype == celltype),]
  plots <- list()

  for (j in 1:nrow(tmpdf)){    
    stages <- as.character(tmpdf[j, 2])
    cluster1 <- as.character(tmpdf[j, 3])
    cluster2 <- as.character(tmpdf[j, 4])
    obToUse <- subset(obtmp, subset = stage == stages & (cluster == cluster1 | cluster == cluster2))
    obToUse <- SetIdent(obToUse, value = obToUse@meta.data$stage) # CHANGING IDENT TO STAGE

    nametoadd <- paste(stages, cluster1, cluster2, sep = "_")
    file = paste("../results/NUCLEAR/", paste(celltype, nametoadd, "cluster zero mtDNA DE genes.xlsx", sep = " "), sep = "")
    if (file.exists(file)){
      DE <- as.data.frame(read.xlsx(file, 1, header=TRUE))
      features <- as.vector(DE[which(DE$p_val <= 0.05 & DE$avg_logFC >= 0.2), 1])      
      #features <- c(features , mtGenes[! mtGenes %chin% features])
      if(length(features) > 1){
        title <- paste(celltype, stages, "Clusters", paste(cluster1, cluster2, sep = "+"), sep = " ")
        plots[[j]] <-  DoHeatmap(obToUse, features = features, size = 4, angle = 90) + 
                          NoLegend() + 
                          ggtitle(title)
      }

      # features <- DE[which(DE$p_val <= 0.05 & DE$avg_logFC >= 0.2), 1]
      # if (length(features) >= 1){
      #   for (k in 1:length(features)){
      #     DE_counter <- DE_counter + 1
      #     tmp <- data.frame(celltype = celltype, clusters = paste(cluster1, cluster2, sep = "_"), feature = features[[k]])
      #     DE_features <- rbind(DE_features, tmp)
      #   }
      # } 
    }
  }
  # filename = paste("Heatmaps", celltype, "mtDNA DE.pdf", sep = " ")
  # pdf(filename, 8, 130)
  # CombinePlots(plots, ncol = 1, legend = "none")
  # dev.off()

  filename = paste("Heatmaps", celltype, "mtDNA DE.pdf", sep = " ")
  pdf(filename, onefile = TRUE)
  for(a in 1:length(plots)){
    do.call("grid.arrange", plots[[a]])
  }
  dev.off()
    
  index <- index + 1
  if(length(celltypes) < index) break
}

print("here is the result:")
DE_features
DE_counter
