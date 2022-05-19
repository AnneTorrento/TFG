library(Seurat)
library(dplyr)
library(data.table)
library(magrittr)
library(tibble)

print("READING NORMALIZED ATLAS")
obj <- readRDS("obj_normalized.rds")

print("READING FEATURES")
features <- readRDS("features_LongPolyA.rds")

print("SUBSETTING OBJECT")
obj = obj[features,]

print("CREATING DATA TABLE")
data <- obj@assays$RNA@data %>% as.matrix %>% t %>% as.data.table

print("CREATING METADATA TABLE")
celltypes <- obj@meta.data$celltype %>% as.data.table
stages <- obj@meta.data$stage %>% as.data.table
md <- cbind(celltypes, stages)
colnames(md) <- c("celltype", "stage")

print("MERGING")
df <- cbind(data, md)
df <- rownames_to_column(df, "sample")

print("SAVING")
write.table(df, "table.txt")

celltypes <- unique(df$celltype)
saveRDS(celltypes, file = "celltypes.rds")
