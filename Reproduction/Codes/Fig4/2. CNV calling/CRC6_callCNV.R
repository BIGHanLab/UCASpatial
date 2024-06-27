# CNV evaluation in CRC6

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
setwd("/data/xy/Spatial_transcriptome/eWEIDE/Reproduction/Data/Fig4_5/clone_region/")
# Create the inferCNV object from the 3 files, without the reference group and mitochondria
NC_CRC6_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="./CRC6_count_matrix.tsv",
                                                 gene_order_file="./siCNV_GeneOrderFile.tsv",
                                                 annotations_file="./Region_CRC6.tsv",
                                                 delim="\t",
                                                 ref_group_names="Remaining",
                                                 chr_exclude = c("chrM","chrX","chrY"))
saveRDS(NC_CRC6_infCNV,file="./CRC6_infCNV_clone_group_2.RDS")

#Run inferCNV
NC_CRC6_infCNV_res_spot = infercnv::run(NC_CRC6_infCNV,
                                        cutoff=0.1, #(see infercnv::run documentation)
                                        out_dir="./SPCNVoutputs_CRC6_clone_group_2",
                                        cluster_by_groups=TRUE, #unsupervised analysis
                                        HMM = TRUE,k_obs_groups = 3,
                                        denoise=TRUE) #denoising applies noise reduction for the plot
saveRDS(NC_CRC6_infCNV_res_spot,file="./CRC6_infCNV_clone_group_2_res_cluster_level_3.RDS")