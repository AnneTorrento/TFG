.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(tidyr)
library(xlsx)

##############################################################################################################
print("READING DATA")

obj <- readRDS("obj_normalized_clusters.rds")
signifClust <- readRDS("signifDiffClusters.rds")
all_celltypes <- readRDS("celltypes.rds")
celltypes <- unique(signifClust$celltype)
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


dode <- function(obtmp,celltype,ident="cluster",nametoadd,namestotest,listtotest){
  # This function writes the excel file of differential expressed genes
  name = celltype
  Idents(obtmp) = ident

  d = FindAllMarkers(obtmp,only.pos=T,min.pct=0.25, logfc.threshold=0,return.thresh=1)
  
  if(dim(d) > 0) { ############################## ADDED IF HERE ##############################

    ptable = getgenes(d,listtotest[[1]])

    start=2
    while(nrow(ptable)<1 && start<length(namestotest)+1){start=start+1;ptable = getgenes(d,listtotest[[start]])}

    if(start<length(namestotest)+1){
      ptable <- ptable %>% drop_na()
      ptable[is.na(ptable)] = 0
      ptable$pathway = namestotest[1]
      ptable$cell_type = celltype
      dataset_name <- paste(stages, cluster1, cluster2, sep = "_")
      ptable$stage_clusters = dataset_name
      for (i in start:length(namestotest)){
        pathway = namestotest[i]
        mygenes = listtotest[[i]]
        ptabletmp = getgenes(d,listtotest[[i]])
        if (length(rownames(ptabletmp))>0){
          ptabletmp <- ptabletmp %>% drop_na()
          ptabletmp[is.na(ptabletmp)] = 0
          ptabletmp$pathway = namestotest[i]
          ptabletmp$cell_type = celltype
          ptabletmp$stage_clusters = dataset_name
          ptable = rbind(ptable,ptabletmp)
        }
      }
      
      list_of_datasets <- list("All_DE_genes" = d,"mtDNA"=getgenes(d,mtGenes), "pathways" = ptable)

      sapply(names(list_of_datasets), 
      function(x) write.xlsx(list_of_datasets[[x]], file = paste(celltype, nametoadd, ident,"zero mtDNA DE genes.xlsx")))

    } 

  } ############################## ADDED IF HERE ##############################
}



##############################################################################################################
print("RUNNING DIFFERENTIAL EXPRESSION ANALYSIS")

nucGenes[[1]] = mtGenes
obj.list <- SplitObject(obj, split.by = "celltype")
index <- 1

for (i in 1:length(all_celltypes)){
  if (all_celltypes[i] != celltypes[index]){
    next
  }

  celltype <- celltypes[index]
  obtmp <- obj.list[[i]]
  tmpdf <- signifClust[which(signifClust$celltype == celltype),]

  for (j in 1:nrow(tmpdf)){    
    stages <- as.character(tmpdf[j, 2])
    cluster1 <- as.character(tmpdf[j, 3])
    cluster2 <- as.character(tmpdf[j, 4])
    obToUse <- subset(obtmp, subset = stage == stages & (cluster == cluster1 | cluster == cluster2))
    
    dode(obToUse, celltype, ident="cluster", nametoadd = paste(stages, cluster1, cluster2, sep = "_"), names(nucGenes), nucGenes)

  }
  index <- index + 1
  if(length(celltypes) < index) break
}
