---
layout: page
title: Example Analysis
description: ~
---
This tutorial is the example analysis with UCASpatial on the human colorectal cancer data from [Qi et al, 2022](https://www.nature.com/articles/s41467-022-29366-6). Before runing the tutorial, make sure that the UCASpatial package is installed. Installation instructions see the [link](https://bighanlab.github.io/UCASpatial/installation.html)

## Required input data
`UCASpatial` requires two types of input data:
- spatial transcriptomics data.  
- single cell RNAseq (scRNA-seq) reference data, along with meta information indicating the cell type information for each cell (clust_vr).  

The example data for runing the tutorial can be downloaded in this [page](https://zenodo.org/records/12634581?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImU0OGM4YjBjLTliYzQtNDc4ZC05Njg3LTczYTdhOTExNTliNiIsImRhdGEiOnt9LCJyYW5kb20iOiI0ZjYyYzU1OGYwNzc4M2E2NDQ5YWFiMjZkZmY2OGE4NCJ9.qcl1aa2yo6mmdkZXdHw2exPmfhyA7M5oPS1ZoKXbqVpBdQbSWkagYIGL2U3Dhj2ffSde0rhe5MFiuIRYuPC_Vw) (need replace)
Here are the details about the required data input illustrated by the example datasets. 

### 1. spatial transcriptomics data, e.g.,
```r
st_vis <- readRDS("./Human_CRC_ST_data.rds")
st_vis
An object of class Seurat 
36601 features across 4248 samples within 1 assay 
Active assay: Spatial (36601 features, 2000 variable features)
 3 layers present: counts, data, scale.data
 2 dimensional reductions calculated: pca, umap
 1 image present: image
```
### 2. single cell RNAseq ((scRNA-seq)) data,  e.g.,
```r
sc_ref <- readRDS('./sc_ref_down.rds')
sc_ref
An object of class Seurat 
45217 features across 2913 samples within 7 assays 
Active assay: RNA (27905 features, 0 variable features)
 2 layers present: counts, data
 6 other assays present: integrated, source, protein, transfered_type, transfered_subtype, SCT
 2 dimensional reductions calculated: pca, umap
```
The cell type information for each cell is stored in the metadata of the Seurat object. The variable name here is `UCASpatial_clus_v7`.

```r
clust_vr <- 'UCASpatial_clus_v7'
table(sc_ref$UCASpatial_clus_v7)

              CD8 Tem               CD8 Tex               CD8 Trm 
                  100                   100                    48 
               CD4 Tn                CD4 Tm                   Tfh 
                  100                   100                    15 
             Th1/Th17                  Treg                    NK 
                   34                   100                   100 
                  ILC                  Mast              Monocyte 
                  100                   100                   100 
           Neutrophil             IGF1+ Mac              FN1+ Mac 
                    2                   100                   100 
            PDPN+ Mac                  cDC1                  cDC2 
                   72                    58                   100 
                 cDC3                   pDC            Fibroblast 
                   58                    31                   100 
          FAP+ Myofib         POSTN+ Myofib         Smooth Muscle 
                  100                   100                   100 
             Pericyte           Endothelium                 Glial 
                  100                   100                   100 
                    B                Plasma         Epithelium c1 
                  100                   100                   100 
LEFTY1+ Epithelium c2 TM4SF1+ Epithelium c3   MUC2+ Epithelium c4 
                  100                   100                   100 
 REG1A+ Epithelium c5 
                   95 
```

## Cell Type Deconvolution
```r
library(UCASpatial)
``` 
### Run UCASpatial deconvolution:
```r
UCASpatial_result <- UCASpatial_deconv(sc_ref =sc_ref,st_vis = st_vis,
                                        clust_vr = clust_vr,downsample_n = 0,
                                        meta.filter = F,cos.filter.threshold=0.05,
                                        ent.filter.threshold = 0.5,
                                        weight.filter.threshold = 0.2)

Load required packages...
Step0    Check the variables............................
Warning: Using current path: "/data/xy/Spatial_transcriptome/eWEIDE/20251201_NC_revision_R1/R3Q8_tutorial" as the 'output_path'...
...........
Load the spatial expression matrix (row counts):st_vis@assays$Spatial@counts
Step0    Preprocess the sc_ref data.....................
Step0.1  Meta.cell filter...............................
Step1    Calculate the markers and add weights..........
Step1.1  Calculate the marker genes.....................
Change the idents of sc_ref into 'sc_ref@meta.data$UCASpatial_clus_v7'...
Calculating cluster CD8 Tem
Calculating cluster CD8 Tex
Calculating cluster CD8 Trm
Calculating cluster CD4 Tn
Calculating cluster CD4 Tm
Calculating cluster Tfh
Calculating cluster Th1/Th17
Calculating cluster Treg
Calculating cluster NK
Calculating cluster ILC
Calculating cluster Mast
Calculating cluster Monocyte
Calculating cluster Neutrophil
Calculating cluster IGF1+ Mac
Calculating cluster FN1+ Mac
Calculating cluster PDPN+ Mac
Calculating cluster cDC1
Calculating cluster cDC2
Calculating cluster cDC3
Calculating cluster pDC
Calculating cluster Fibroblast
Calculating cluster FAP+ Myofib
Calculating cluster POSTN+ Myofib
Calculating cluster Smooth Muscle
Calculating cluster Pericyte
Calculating cluster Endothelium
Calculating cluster Glial
Calculating cluster B
Calculating cluster Plasma
Calculating cluster Epithelium c1
Calculating cluster LEFTY1+ Epithelium c2
Calculating cluster TM4SF1+ Epithelium c3
Calculating cluster MUC2+ Epithelium c4
Calculating cluster REG1A+ Epithelium c5
Auto save the marker genes under the path:
/data/xy/Spatial_transcriptome/eWEIDE/20251201_NC_revision_R1/R3Q8_tutorial/cluster_markers.rds
Step1.2  Calculate the entropy-based weight.............
  |=====================================================================| 100%

Warning: Setting row names on a tibble is deprecated.

Step1.3  Calculate the cosine-based weight..............
Warning: x or y has vectors with all zero; consider setting use_nan = TRUE to set these values to NaN or use_nan = FALSE to suppress this warning
Step1.4  Filter the markers.............................
Save the markers........................................
Warning: Layer counts isn't present in the assay object; returning NULL
Step2  Train the nsNMF model............................
[1] "Preparing Gene set"
Warning: Layer counts isn't present in the assay object; returning NULL
Warning in asMethod(object) :
  sparse->dense coercion: allocating vector of size 1.2 GiB
Warning: Layer counts isn't present in the assay object; returning NULL
Normalize the sc_ref matrix...
Initialize the NMF matrices...
[1] "NMF Training..."
[1] "Time to initialize and train NMF model was 0.05mins"
Step3  Calculate the cluster-topic profile..............
Step4  Weighted-NNLS to implement the deconvolution.....
[1] "Deconvoluting spots"
  |=====================================================================| 100%
```

### Quality Control Plot

```r
p <- dot_plot_profiles_fun(
    UCASpatial_result[[1]][[1]]@h,
    UCASpatial_result[[1]][[2]])[2]
print(p)
```
Here is an example output: 
![Example_Dotplot](QC_plot_UCASpatial.png)

### Visualization of cell type distributions
```r
# Extract results
decon_matr <- as.matrix(UCASpatial_result[[2]])
cell_proportions <- decon_matr[,1:(ncol(decon_matr)-1)]
cell_proportions <- cell_proportions/rowSums(cell_proportions)
rownames(cell_proportions) <- colnames(st_vis)

# Prepare data for visualization
decon_df <- cell_proportions %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")

# Add to Seurat object
st_vis@meta.data <- st_vis@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")

# Create assay for cell type proportions
Annotation_assay <- CreateAssayObject(t(cell_proportions))
st_vis@assays$Annotation <- Annotation_assay
st_vis@assays$Annotation@key <- "annatation_"
DefaultAssay(st_vis) <- "Annotation"
# Plot cell type distributions
cell_types_all <- colnames(cell_proportions)
plot_list <- Seurat::SpatialFeaturePlot(
    object = st_vis,
    features = c('FAP..Myofib','CD4.Tn','B'),
    stroke = NA,
    pt.size.factor = 1000,
    alpha = c(0.3, 1),
    min.cutoff = 0.03
)
plot_list
```
Here is an example output: 
![Example_Spatial plot](Cell_type_distributions_UCASpatial.png)
