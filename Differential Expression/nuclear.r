library(Seurat)
library(pvclust)


##############################################################################################################
print("READING DATA")

# obj <- readRDS("obj_normalized.rds")
obj <- readRDS("obj_normalized_clusters.rds")
celltypes <- readRDS("celltypes.rds")
mtGenes <- readRDS("features_LongPolyA.rds")
nucGenes <- readRDS("genesets_symb.rds")  #List of orthologous genes from the nucleus that are imported into the mitochondria. 
                                          #The names are the pathways' names.


##############################################################################################################
print("DEFINING FUNCTIONS")

getgenes<- function(DE, mygenes){
  # This function subsets a set of genes from the original differential expression list
  # correct it for multiple testing
  de = DE[rownames(DE) %in% mygenes,]
  if (length(rownames(de))<1){return(de)}
  de$p_val_adj_FDR = p.adjust(de$p_val, method = "fdr") 
  de$p_val_adj_BH = p.adjust(de$p_val, method = "BH") 
  return(de)  
}


dode <- function(obtmp,celltype,ident="severity",nametoadd,namestotest,listtotest){
  # This function writes the excel file of differential expressed genes
  name = celltype
  Idents(obtmp) = ident
  d = FindAllMarkers(obtmp,only.pos=T,min.pct=0.25, logfc.threshold=0,return.thresh=1)
  ptable = getgenes(d,listtotest[[1]])
  start=2
  while(nrow(ptable)<1 && start<length(namestotest)+1){start=start+1;ptable = getgenes(d,listtotest[[start]])}
  if(start<length(namestotest)+1){
    ptable <- ptable %>% drop_na()
    ptable[is.na(ptable)] = 0
    ptable$pathway = namestotest[1]
    ptable$cell_type = celltype
    ptable$dataset = dataset_name
    for (i in start:length(namestotest)){
      pathway = namestotest[i]
      mygenes = listtotest[[i]]
      ptabletmp = getgenes(d,listtotest[[i]])
      if (length(rownames(ptabletmp))>0){
        ptabletmp <- ptabletmp %>% drop_na()
        ptabletmp[is.na(ptabletmp)] = 0
        ptabletmp$pathway = namestotest[i]
        ptabletmp$cell_type = celltype
        ptabletmp$dataset = dataset_name
        ptable = rbind(ptable,ptabletmp)
      }
    }
    
    list_of_datasets <- list("All_DE_genes" = d,"mtDNA"=getgenes(d,mito.genes), "pathways" = ptable)
    write.xlsx(list_of_datasets, file = paste(nametoadd,celltype,ident,"zero mtDNA DE genes.xlsx"))
  }
}



##############################################################################################################
print("ADDING HIERARCHICAL (PVCLUST) CLUSTERS TO OBJECT METADATA")
nClusters <- c(5, 6, 4, 4, 5,
               4, 5, 5, 4, 5,
               5, 5, 7, 5, 4,
               5, 7, 4, 7, 5,
               4, 5, 8, 5, 4,
               5, 4, 5, 7, 4,
               5, 4, 4, 2, 4, 
               2, 2)

cluster <- data.frame()

for (i in 1:length(celltypes)){
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- as.data.frame(cutree(result$hclust, k=nClusters[i]))
  rownames(tmp) <- NULL
  cluster <- rbind(cluster, tmp)
}

colnames(cluster) <- "cluster"

obj@meta.data <- cbind(obj@meta.data, cluster)
names(obj@meta.data[ncol(obj@meta.data)]) <- colnames(cluster) #or <- "cluster"

saveRDS(obj, file = "obj_normalized_clusters.rds")



##############################################################################################################
print("CREATING DATAFRAME WITH CELLTYPES AND STAGES WHERE TWO CELL CLUSTERS DIFFER BY EXPRESSION") #(AFTER WILCOXON WITH HIERARCHICAL CLUSTERING)
stage_cluster <- data.frame(matrix(ncol = 4, nrow = 0))

for(i in 1:length(celltypes)){
  print(paste("LOADING DATA", celltypes[i], sep = " "))
  result <- readRDS(paste("../results/HIERARCHICAL/PVCLUST/pvclust_", celltypes[i],".rds", sep = ""))
  tmp <- my_data[which(my_data$celltype == celltypes[i]),]
  tmp$cluster <- cutree(result$hclust, k=nClusters[i])
  tmp$AvgExprs <- rowMeans(tmp[,2:11])
  
  for(j in 1:length(unique(tmp$stage))){
    stage <- which(tmp$stage == unique(tmp$stage)[j])
    
    if(stage){
      tmp_stage <- as.data.frame(cbind(tmp[stage,15],
                                       tmp[stage,13],
                                       tmp[stage,14]))
      colnames(tmp_stage) <- c("AvgExprs", "stage", "cluster")
      tmp_stage <- transform(tmp_stage, AvgExprs = as.numeric(AvgExprs))
      
      if(length(unique(tmp_stage$cluster)) >= 2){
        test <- tmp_stage %>% wilcox_test(AvgExprs ~ cluster)
      }
      
      if(length(unique(tmp_stage$cluster)) == 2){
        if (test$p <= 0.01){
          newrow <- c(as.character(celltypes[i]), unique(tmp_stage$stage), test$group1, test$group2)
          stage_cluster <- rbind(stage_cluster, newrow)
        }
      }
      
      else if(length(unique(tmp_stage$cluster)) > 2){
        if (test$p.adj <= 0.01){
          for (k in 1:length(test$group1)){
            newrow <- c(as.character(celltypes[i]), unique(tmp_stage$stage), test$group1[k], test$group2[k])
            stage_cluster <- rbind(stage_cluster, newrow)
          }
        }
      }
      
    }
  }
}

colnames(stage_cluster) <- c("celltype", "stage", "cluster 1", "cluster 2")
saveRDS(stage_cluster, file = "signifDiffClusters_perStage.rds")


##############################################################################################################
print("RUNNING DIFFERENTIAL EXPRESSION ANALYSIS")

nucGenes[[1]] = mtGenes
obj.list <- SplitObject(obj, split.by = "celltype")

stagelist <- 

clusterlist <-

for (i in 1:length(celltypes)){
  celltype <- celltypes[i]
  obtmp <- obj.list[[i]]
  stages <- unique(obtmp@meta.data$stage)

  for (j in 1:length(stages)){
    obtmp_stage.list <- SplitObject(obtmp, split.by = "stage")
    obtmp_stage <- obtmp_stage.list[[j]]
    clusters <- unique(obtmp_stage@meta.data$cluster)
    
    for (k in 1:length()){
      dode(obtmp, celltype, ident="severity", "Nuc", names(nucGenes), nucGenes)
    }
  }
}

