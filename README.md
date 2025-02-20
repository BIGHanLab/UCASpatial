# UCASpatial
## Ultra-resolution Deconvolution of Spatial Transcriptomics Unveils Spatiotemporal Cellular Dynamics in Complex Microenvironments
![image](https://github.com/BIGHanLab/UCASpatial/assets/167292686/96adbd00-41cd-49ee-a0d8-2e0fc52d33ba)

Here, we introduce an ultra-resolution ST deconvolution algorithm (UCASpatial) that improves the resolution of mapping cell subpopulations to spatial locations by leveraging the contribution of genes indicative of cell identity through entropy-based weighting. Using both in silico and real ST datasets, we demonstrate that UCASpatial improves the robustness and accuracy in identifying low-abundant cell subpopulations and distinguishing transcriptionally heterogeneous cell subpopulations. 

## Installation (R)
You can install the released version of UCASpatial from GitHub by:
```R
devtools::install_github('https://github.com/BIGHanLab/UCASpatial/')
```
You can also install the released version of UCASpatial from GitHub by download the 'UCASpatial_v1.R', and then source it in Rstudio on your own working direction.
```R
source('/DataPath/UCASpatial_v1.R')
```
## Dependencies
* R >= 3.6.1
* Seurat >= 3.1.4
* RcppML >= 0.5.6
## How to use UCASpatial
See 'Reproduction' for the reproduction of the main figures of the paper.
To briefly use UCASpatial, first need to load the spatial transcriptomics data and the referenced scRNA/snRNA-seq data:
```R
st_vis <- readRDS('/datapath/ST_data.rds')
sc_ref <- readRDS('/datapath/SC_data.rds') # cluster info in 'sc.ref$clust_vr'
```
Run UCASpatial to map cell subpopulations to spatial locations
```R
UCASpatial_result <- UCASpatial_deconv(
      sc_ref = sc_ref,
      st_vis = st_vis,
      clust_vr = clust_vr)
```
UCASpatial also support the analysis on 10X Visium HD datasets:
```R
# Run UCASpatial for 10X HD data
UCASpatial_result <- UCASpatial_HD_deconv(
      sc_ref = sc_ref,
      st_vis = st_vis,
      spatial.assay='Spatial.008um',
      clust_vr = clust_vr)
```

```R
# Evaluation
p <- dot_plot_profiles_fun(UCASpatial_result[[1]][[1]]@h,UCASpatial_result[[1]][[2]])[2]
p
```
```R
# Visualization
decon_matr <- as.matrix(UCASpatial_result[[2]])[,1:(ncol(UCASpatial_result[[2]])-1)]
decon_matr <- decon_matr/rowSums(decon_matr)
rownames(decon_matr) <- colnames(st_vis)
decon_df <- decon_matr %>%
      data.frame() %>%
      tibble::rownames_to_column("barcodes")
st_vis@meta.data <- st_vis@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(decon_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
Annotation_assay <- CreateAssayObject(t(decon_matr))
st_vis@assays$Annotation <- Annotation_assay
cell_types_all <- colnames(decon_matr)
Seurat::SpatialFeaturePlot(
      object = st_vis,
      features = cell_types_all,stroke = NA,alpha = c(0.3,1),
      min.cutoff = 0.03)
```
## UCASpatial for python
The core functionality of UCASpatial is already available in Python.
You can utilize the released version of UCASpatial for Python from GitHub by downloading 'UCASpatial_python_v1.0.py' and importing it. Application examples can be found under the Python path of the current GitHub repository.

For using UCASpatial in Python:
```python
import UCASpatial_python_v1.0

# Load the data
sc_ref = sc.read("your sc_ref path/sc.h5ad")
st_vis = sc.read("your st_vis path/st.h5ad")

X_norm = sc.pp.normalize_total(sc_ref,target_sum=1,inplace=False)['X']
sc_ref.layers['data'] = X_norm

# Initialization and Execution
ucas = UCASpatial_ds_R1.UCASpatial(
    sc_ref=sc_ref,
    st_vis=st_vis,
    clust_vr='your clust_vr',
    meta_filter=False,
    random_seed=12345
)
result = ucas.run()

# Get the results
nmf_components = result['nmf_components']
cell_proportions = result['proportions']
cell_proportions.to_csv('cell_proportions.csv')
```

## Issues
All feedback, bug reports, and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R!

## How to cite UCASpatial
Xu, Y. et al. Ultra-resolution Deconvolution of Spatial Transcriptomics Unveils Spatiotemporal Cellular Dynamics in Complex Microenvironments. bioRxiv, 2024.2007.2005.602200 (2024).


