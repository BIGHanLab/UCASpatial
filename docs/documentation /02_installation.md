---
layout: page
title: Installation
description: ~
---

`UCASpatial` is implemented as an R package, which can be installed from GitHub by:

### Dependencies
* R >= 3.6.1
* Seurat >= 3.1.4
* RcppML >= 0.5.6

#### 1. Install `devtools` if necessary
```r
install.packages('devtools')
```

#### 2. Install `UCASpatial` 
You can install the released version of UCASpatial from GitHub by:
```R
devtools::install_github('https://github.com/BIGHanLab/UCASpatial/tree/UCASpatial_v1.4.0')
```
You can also install the released version of UCASpatial from GitHub by download the 'UCASpatial_v1.R', and then source it in Rstudio on your own working direction.
```R
source('/DataPath/UCASpatial_v1.R')
```

#### 3. Load package
```r
library(UCASpatial)
```
