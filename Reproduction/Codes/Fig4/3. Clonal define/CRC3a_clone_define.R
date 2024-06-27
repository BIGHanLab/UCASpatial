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
CRC3a_st<-readRDS("CRC3a_region.RDS")
#Use read.dendrogram() to import the dendogram file
CRC3a_st_for_clustering <- read.dendrogram(file = "CRC3a_infercnv.21_denoised.observations_dendrogram.txt")

#Convert to a phylo object using as.phylo()
CRC3a_st_for_clustering_phylo <- as.phylo(CRC3a_st_for_clustering)

#Use subtrees() to enable further interaction with the dendrogram
my.subtrees = subtrees(CRC3a_st_for_clustering_phylo)  # subtrees() to subset

#Output an image to visualize all of the dengdrogram nodes 
png("CRC3a_st_forclustering_phylo.png",width=10000,height=2500, res = 300)
plot(CRC3a_st_for_clustering_phylo,show.tip.label = FALSE)
nodelabels(text=1:CRC3a_st_for_clustering_phylo$Nnode,node=1:CRC3a_st_for_clustering_phylo$Nnode+Ntip(CRC3a_st_for_clustering_phylo))
dev.off()



#Selecting Node 1
Node1 <- SelectingSubTreeData(my.subtrees, 337)

#Selecting Node 2
Node2 <- SelectingSubTreeData(my.subtrees, 2)

#Selecting Node 3
#Node3 <- SelectingSubTreeData(my.subtrees, 2)

#Merging All Nodes Together
MergedNodes <- rbind(Node1, Node2)
#MergedNodes <- Node1

#Note, since the "CRC3aitional spot" for Clone 3 did not have a node, identify it by joining the original annotations to the MergedNodes dataframe
#MergedNodes <- full_join(MergedNodes, Celltype_Anno_CRC3a_st)


#Then drop the Histology column
#MergedNodes <- MergedNodes[,1:2]

#This then renames "Node" to "Clone"
MergedNodes$Node <- ifelse(MergedNodes$Node == "Node_337", "Clone1",
                           ifelse(MergedNodes$Node == "Node_2", "Clone2",MergedNodes$Node))

#Show the first part of the finalized dataframe for this step
head(MergedNodes)


#Copy the MergedNodes dataframe to a new dataframe called "ForLoupeBrowser"
ForLoupeBrowser <- MergedNodes

#Replace the "." in the current barcode dataframe to "-" (required for LoupeBrowser)
#ForLoupeBrowser$Barcode <- gsub("\\.", "\\-", ForLoupeBrowser$Barcode)
#ForLoupeBrowser$Barcode <- trimws(substr(ForLoupeBrowser$Barcode, 8, 100))
#save Loupe object
write.csv(ForLoupeBrowser, "CRC3a_st_Clones_ForLoupeBrowser.csv", row.names = FALSE)
#CRC3a_st@meta.data<-CRC3a_st@meta.data
CRC3a_st@meta.data$Barcode<-rownames(CRC3a_st@meta.data)
CRC3a_st@meta.data<-merge(CRC3a_st@meta.data,ForLoupeBrowser,by="Barcode",all = T)
rownames(CRC3a_st@meta.data) <- CRC3a_st@meta.data$Barcode
CRC3a_st@meta.data<-CRC3a_st@meta.data[colnames(CRC3a_st),]
rownames(CRC3a_st@meta.data) <- CRC3a_st@meta.data$Barcode
SpatialPlot(CRC3a_st,group.by = c("Node"))
library(Seurat)
#CRC3a_st$Node.y[CRC3a_st$Node.y %in% c("Clone4","Clone5")] <- "Clone4"
#CRC3a_st$Region[CRC3a_st$Node %in% c("Clone1","Clone2")] <- "Clone1"
CRC3a_st$Region[CRC3a_st$Node %in% c("Clone1")] <- "Clone1"
CRC3a_st$Region[CRC3a_st$Node %in% c("Clone2")] <- "Normal_Epithelial"
#CRC3a_st$Region[CRC3a_st$Node %in% c("Clone4")] <- "Clone3"
pdf("clone_distribution_CRC3a_st.pdf",width = 6,height = 5)
SpatialPlot(CRC3a_st,group.by = c("Region"),cols = c("#F6B272","#B2DF8A","white"),image.alpha = 0)
dev.off()
saveRDS(CRC3a_st,file = "../clone_region/CRC3a_st_clone.RDS")
