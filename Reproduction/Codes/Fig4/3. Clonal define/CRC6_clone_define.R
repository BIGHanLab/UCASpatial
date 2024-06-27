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
CRC6_st<-readRDS("CRC6_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC6_st_for_clustering <- read.dendrogram(file = "CRC6_infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC6_st_for_clustering_phylo <- as.phylo(CRC6_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC6_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("CRC6_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC6_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC6_st_for_clustering_phylo$Nnode,node=1:CRC6_st_for_clustering_phylo$Nnode+Ntip(CRC6_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 1989)

#Selecting Node 2
Node2 <- SelectingSubTreeData(my.subtrees, 738)

#Selecting Node 3
Node3 <- SelectingSubTreeData(my.subtrees, 2)

#Selecting Node 3
#Node4 <- SelectingSubTreeData(my.subtrees, 2)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2,Node3)
#MergedNodes <- Node1

#Note, since the "CRC6itional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC6_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_738", "Clone1",
                           ifelse(MergedNodes$Node == "Node_1989", "Clone2",
                                  ifelse(MergedNodes$Node == "Node_2", "Normal_Epithelial",MergedNodes$Node)))

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC6_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC6_st@meta.data<-CRC6_st@meta.data
CRC6_st@meta.data$Barcode<-rownames(CRC6_st@meta.data)
CRC6_st@meta.data<-merge(CRC6_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC6_st@meta.data) <- CRC6_st@meta.data$Barcode
CRC6_st@meta.data<-CRC6_st@meta.data[colnames(CRC6_st),]
rownames(CRC6_st@meta.data) <- CRC6_st@meta.data$Barcode
SpatialPlot(CRC6_st,group.by = c("Node"))
library(Seurat)
#CRC6_st$Node.y[CRC6_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC6_st$Region[CRC6_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC6_st$Region[CRC6_st$Node %in% c("Clone1")] <- "Clone1"
CRC6_st$Region[CRC6_st$Node %in% c("Clone2")] <- "Remaining"
CRC6_st$Region[CRC6_st$Node %in% c("Normal_Epithelial")] <- "Normal_Epithelial"
#CRC6_st$Region[CRC6_st$Node %in% c("Clone4")] <- "Clone3"
pdf("clone_distribution_CRC6_st.pdf",width = 6,height = 5)
SpatialPlot(CRC6_st,group.by = c("Region"),cols = c("#F6B272","#B2DF8A","white"),image.alpha = 0)
dev.off()
saveRDS(CRC6_st,file = "../clone_region/CRC6_st_clone.RDS")
