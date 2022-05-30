library(xlsx)
library(gplots)
library(RColorBrewer)

ectoderm <- read.xlsx("heatmap_data.xlsx", 1, header = TRUE)
endoderm <- read.xlsx("heatmap_data.xlsx", 2, header = TRUE)
mesoderm <- read.xlsx("heatmap_data.xlsx", 3, header = TRUE)

ecto <- ectoderm[,-1]
endo <- endoderm[,-1]
meso <- mesoderm[,-1]
rownames(ecto) <- ectoderm[,1]
rownames(endo) <- endoderm[,1]
rownames(meso) <- mesoderm[,1]

m_ecto <- as.matrix(ecto)
m_endo <- as.matrix(endo)
m_meso <- as.matrix(meso)

colors <- c("#F2F2F2", "#ADD8E6", "#1AA7EC")

png("heat_ecdo.png")
heatmap.2(m_ecto, dendrogram = "none", 
          Rowv = FALSE, Colv = FALSE, trace = "none", 
          breaks=seq(-2,1,1), col=colors,
          margins = c(10, 14))
dev.off()

png("heat_endo.png")
heatmap.2(m_endo, dendrogram = "none", 
          Rowv = FALSE, Colv = FALSE, trace = "none", 
          breaks=seq(-2,1,1), col=colors,
          margins = c(10, 14))
dev.off()

png("heat_meso.png")
heatmap.2(m_meso, dendrogram = "none", 
          Rowv = FALSE, Colv = FALSE, trace = "none", 
          breaks=seq(-2,1,1), col=colors,
          margins = c(10, 14))
dev.off()
