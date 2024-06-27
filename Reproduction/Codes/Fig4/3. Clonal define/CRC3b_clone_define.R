# 3. clonal region annotation
setwd("/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Data/Fig4_5/clonal_define/")
library("devtools")
library("BiocManager")
library("infercnv")
library("tidyverse")
library("Seurat")
library("phylogram")
library("ape")
library("hdf5r")
library("SpatialInferCNV")
library(ape)
library(pacman)
p_unload(Seurat)
p_load(Seurat)
CRC3b_st<-readRDS("CRC3b_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC3b_st_for_clustering <- read.dendrogram(file = "CRC3b_infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC3b_st_for_clustering_phylo <- as.phylo(CRC3b_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC3b_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("CRC3b_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC3b_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC3b_st_for_clustering_phylo$Nnode,node=1:CRC3b_st_for_clustering_phylo$Nnode+Ntip(CRC3b_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 1)

#Selecting Node 2
#Node2 <- SelectingSubTreeData(my.subtrees, 2)

#Selecting Node 3
#Node3 <- SelectingSubTreeData(my.subtrees, 2)

#Merging All Nodes Together
MergedNodes <- Node1
#MergedNodes <- Node1

#Note, since the "CRC3bitional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC3b_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_1", "Clone1",MergedNodes$Node)

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC3b_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC3b_st@meta.data<-CRC3b_st@meta.data
CRC3b_st@meta.data$Barcode<-rownames(CRC3b_st@meta.data)
CRC3b_st@meta.data<-merge(CRC3b_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC3b_st@meta.data) <- CRC3b_st@meta.data$Barcode
CRC3b_st@meta.data<-CRC3b_st@meta.data[colnames(CRC3b_st),]
rownames(CRC3b_st@meta.data) <- CRC3b_st@meta.data$Barcode
SpatialPlot(CRC3b_st,group.by = c("Node"))
library(Seurat)
#CRC3b_st$Node.y[CRC3b_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC3b_st$Region[CRC3b_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC3b_st$Region[CRC3b_st$Node %in% c("Clone1")] <- "Clone1"
#CRC3b_st$Region[CRC3b_st$Node %in% c("Clone2")] <- "Clone2"
#CRC3b_st$Region[CRC3b_st$Node %in% c("Clone4")] <- "Clone3"
pdf("clone_distribution_CRC3b_st.pdf",width = 6,height = 5)
SpatialPlot(CRC3b_st,group.by = c("Region"),cols = c("#F6B272","#96C2D6","white"),image.alpha = 0)
dev.off()
saveRDS(CRC3b_st,file = "../clone_region/CRC3b_st_clone.RDS")
