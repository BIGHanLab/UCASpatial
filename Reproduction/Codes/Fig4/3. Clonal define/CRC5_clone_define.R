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
CRC5_st<-readRDS("CRC5_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC5_st_for_clustering <- read.dendrogram(file = "CRC5_infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC5_st_for_clustering_phylo <- as.phylo(CRC5_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC5_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("CRC5_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC5_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC5_st_for_clustering_phylo$Nnode,node=1:CRC5_st_for_clustering_phylo$Nnode+Ntip(CRC5_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 1493)

#Selecting Node 2
Node2 <- SelectingSubTreeData(my.subtrees, 1060)

#Selecting Node 3
Node3 <- SelectingSubTreeData(my.subtrees, 517)
Node5 <- SelectingSubTreeData(my.subtrees, 746)

#Selecting Node 3
Node4 <- SelectingSubTreeData(my.subtrees, 2)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2,Node3,Node4,Node5)
#MergedNodes <- Node1

#Note, since the "CRC5itional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC5_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_1493", "Clone1",
                           ifelse(MergedNodes$Node == "Node_1060", "Clone2",
                                  ifelse(MergedNodes$Node == "Node_517", "Clone3",
                                         ifelse(MergedNodes$Node == "Node_2", "Normal_Epithelial",
                                                ifelse(MergedNodes$Node == "Node_746", "Clone4",MergedNodes$Node)))))

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC5_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC5_st@meta.data<-CRC5_st@meta.data
CRC5_st@meta.data$Barcode<-rownames(CRC5_st@meta.data)
CRC5_st@meta.data<-merge(CRC5_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC5_st@meta.data) <- CRC5_st@meta.data$Barcode
CRC5_st@meta.data<-CRC5_st@meta.data[colnames(CRC5_st),]
rownames(CRC5_st@meta.data) <- CRC5_st@meta.data$Barcode
SpatialPlot(CRC5_st,group.by = c("Node"))
library(Seurat)
#CRC5_st$Node.y[CRC5_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC5_st$Region[CRC5_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC5_st$Region[CRC5_st$Node %in% c("Clone1")] <- "Remaining"
CRC5_st$Region[CRC5_st$Node %in% c("Clone2")] <- "Clone1"
CRC5_st$Region[CRC5_st$Node %in% c("Clone3")] <- "Clone2"
CRC5_st$Region[CRC5_st$Node %in% c("Clone4")] <- "Clone3"
CRC5_st$Region[CRC5_st$Node %in% c("Normal_Epithelial")] <- "Normal_Epithelial"
#CRC5_st$Region[CRC5_st$Node %in% c("Clone4")] <- "Clone3"
pdf("clone_distribution_CRC5_st.pdf",width = 6,height = 5)
SpatialPlot(CRC5_st,group.by = c("Region"),cols = c("#F6B272","#96C2D6","#AA2E2E","#B2DF8A","white"),image.alpha = 0)
dev.off()
saveRDS(CRC5_st,file = "../clone_region/CRC5_st_clone.RDS")
