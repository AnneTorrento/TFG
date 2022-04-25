.libPaths("/gpfs0/biores/users/mishmarlab/Anne/R-3.6.3/library")
library(Seurat)
library(ggplot2)
library(reshape)


####################################################################################
#DEFINITIONS
####################################################################################
filtermitozero <- function(ob,n=0,genes=mito.genes){
  # Filter seurat object - exclude cells with more than n zeros in mtdna genes.
  # input - ob = seurat object, mito.genes= vector of mtDNA genes.
  # output - filtered seurat object.
  df<- as.data.frame(ob@assays$RNA@data[genes,])
  exclude = names(df)[colSums(pmax(df == 0)) > n]
  ob = subset(ob, cells = colnames(ob)[!colnames(ob) %in% exclude])
  return(ob)
}

#Protein-coding mitochondrial genes
features <- c("mt-Cytb", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd4l", "mt-Nd5", "mt-Nd6", "mt-Atp6", "mt-Atp8", "mt-Co1", "mt-Co2", "mt-Co3")

#Protein-coding mitochondrial genes with long poly-A tail
features_LongPolyA <- c("mt-Cytb", "mt-Nd1", "mt-Nd2", "mt-Nd3", "mt-Nd4", "mt-Nd5", "mt-Atp6", "mt-Co1", "mt-Co2", "mt-Co3")


####################################################################################
print("PART 1: READING RAW DATA")

data_dir <- "../E-MTAB-6967.processed.1/atlas_data/atlas"
list.files(data_dir)

atlas.data <- Read10X(data.dir = data_dir)
atlas <- CreateSeuratObject(counts = atlas.data, project = "bgu", min.cells = 3, min.features = 200)

atlas[["percent.mt"]] <- PercentageFeatureSet(atlas, pattern = "^mt-")

save(atlas.data, atlas, file = "1atlas_raw.RData")

####################################################################################
print("PART 2: ADDING METADATA")
meta <- read.csv("../E-MTAB-6967.processed.1/atlas_data/atlas/meta.csv")
rownames(meta) <- meta[,1]
meta[,1] <- NULL
atlas@meta.data <-cbind(atlas@meta.data, meta)
names(atlas@meta.data[5:ncol(atlas@meta.data)]) <- colnames(meta)

saveRDS(atlas, file = "ATLAS2_raw_metadata.rds")


####################################################################################
print("PART 3: FILTERING MITOCHONDRIAL ZEROS")
atlas <- filtermitozero(atlas, 0, features_LongPolyA)
saveRDS(atlas, file = "ATLAS3_filtered.rds")


####################################################################################
print("PART 4.1: REMOVING MIXED GASTRULATION CELLS")
atlas <- subset(atlas, subset = stage != "mixed_gastrulation")

print("PART 4.2: REMOVING <NA> CELLTYPE CELLS")
atlas <- subset(atlas, subset = celltype != "<NA>")

saveRDS(atlas, file = "ATLAS4.rds")


####################################################################################
print("PART 5: NORMALIZING DATA")
atlas <- NormalizeData(atlas)
saveRDS(atlas, file = "ATLAS5_normalized.rds")


####################################################################################
print("PART 6: SCALING DATA")
atlas <- ScaleData(atlas)
saveRDS(atlas, file = "ATLAS6_scaled.rds")


####################################################################################
print("PART 7: CREATING MODIFIABLE OBJECTS")
Raw <- atlas@assays$RNA@counts
obj <- CreateSeuratObject(counts = Raw, project = "bgu")
obj@meta.data$celltype <- atlas@meta.data$celltype
obj@meta.data$stage <- atlas@meta.data$stage

obj <- NormalizeData(obj)
saveRDS(obj, file = "obj_normalized.rds")

obj <- ScaleData(obj)
saveRDS(obj, file = "obj_scaled.rds")


