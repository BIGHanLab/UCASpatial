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
CRC4_st<-readRDS("CRC4_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC4_st_for_clustering <- read.dendrogram(file = "CRC4_infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC4_st_for_clustering_phylo <- as.phylo(CRC4_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC4_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("CRC4_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC4_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC4_st_for_clustering_phylo$Nnode,node=1:CRC4_st_for_clustering_phylo$Nnode+Ntip(CRC4_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 1988)
Node5 <- SelectingSubTreeData(my.subtrees, 993)

#Selecting Node 2
Node2 <- SelectingSubTreeData(my.subtrees, 816)

#Selecting Node 3
Node3 <- SelectingSubTreeData(my.subtrees, 2)

#Selecting Node 3
Node4 <- SelectingSubTreeData(my.subtrees, 2086)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2,Node3,Node4,Node5)
#MergedNodes <- Node1

#Note, since the "CRC4itional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC4_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_993", "Clone1",
                           ifelse(MergedNodes$Node == "Node_816", "Clone2",
                                  ifelse(MergedNodes$Node == "Node_2", "Clone3",
                                         ifelse(MergedNodes$Node == "Node_2086", "Remaining",
                                                ifelse(MergedNodes$Node == "Node_1988", "Clone1",MergedNodes$Node)))))

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC4_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC4_st@meta.data<-CRC4_st@meta.data
CRC4_st@meta.data$Barcode<-rownames(CRC4_st@meta.data)
CRC4_st@meta.data<-merge(CRC4_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC4_st@meta.data) <- CRC4_st@meta.data$Barcode
CRC4_st@meta.data<-CRC4_st@meta.data[colnames(CRC4_st),]
rownames(CRC4_st@meta.data) <- CRC4_st@meta.data$Barcode
SpatialPlot(CRC4_st,group.by = c("Node"))
library(Seurat)
#CRC4_st$Node.y[CRC4_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC4_st$Region[CRC4_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC4_st$Region[CRC4_st$Node %in% c("Clone1")] <- "Clone1"
CRC4_st$Region[CRC4_st$Node %in% c("Clone2")] <- "Clone2"
CRC4_st$Region[CRC4_st$Node %in% c("Clone3")] <- "Clone3"
CRC4_st$Region[CRC4_st$Node %in% c("Remaining")] <- "Remaining"
#CRC4_st$Region[CRC4_st$Node %in% c("Clone4")] <- "Clone3"
pdf("clone_distribution_CRC4_st.pdf",width = 6,height = 5)
SpatialPlot(CRC4_st,group.by = c("Region"),cols = c("#F6B272","#96C2D6","#AA2E2E","white"),image.alpha = 0)
dev.off()
saveRDS(CRC4_st,file = "../clone_region/CRC4_st_clone.RDS")
