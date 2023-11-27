# Cardinal Scripts

A small collection of scripts written to streamline data processing/analysis in Cardinal.

## Installation
#### Install from Github
```
devtools::install_github("gtluu/cardinalscripts")
```

#### Load from Local Folder
```
devtools::load_all("/path/to/cardinalscripts")
```

## Dependencies
- R (>= 4.3.1)
- Cardinal (>= 3.4.1)
- stringr(>= 1.5.0)
- fuzzyjoin (>= 0.1.6)
- dplyr (>= 1.1.3)
- ggplot2 (>= 3.4.4)

## Usage

### updateMetadata

Use this function to add region of interest (ROI) and/or condition/treatment information based on x-y coordinates. An example Spot List (exported from Bruker flexImaging) and ROI table can be found in the ```examples``` folder.

```
data <- readMSIData('data.imzML')
data <- updateMetadata(data, spotsFile="spots.txt", roiFile="rois.csv")
```

### subsetData

Use this function to remove any unwanted ROIs from your dataset. This will remove any pixel from that ROI from the ```MSImagingExperiment``` object.

```
data <- readMSIData('data.imzML')  # dataset with 4 ROIs
roi1 <- subset(data, c('roi2', 'roi3', 'roi4'))
roi2 <- subset(data, c('roi1', 'roi3', 'roi4'))
roi3 <- subset(data, c('roi1', 'roi2', 'roi4'))
roi4 <- subset(data, c('roi1', 'roi2', 'roi3'))
```

### ionImageReport

Use this function to generate a PDF containing ion images for each feature from a vector of features.

```
ionImageReport(data, mz=mz(data), filename='ion_images')
```

### sscImageReport

Use this function to generate a PDF containing segmentation maps for each model that was used for multivariate segmentation via spatial shrunken centroids (i.e. each combination of r, k, and s parameters).

```
sscImageReport(ssc, filename='segmentation_maps')
```

### optimizeSSCParams

Use this function to assist in optimize the spatial neighborhood radius (r), initial number of segments (k), and sparsity (s) parameters used for running spatial shrunken centroids segmentation.

```
ssc <- spatialShrunkenCentroids(preprocessed_data, r=c(1, 2), k=c(5, 10, 20), s=c(0, 3, 6, 9))
ssc_params <- optimizeSSCParams(ssc)
```

### getStatisticTable

Use this function to get a ```data.frame``` of values for a given set of r, k, and s parameters that can be used to tell which features contribute significantly to specific classes in a given spatial shrunken centroids segmentation model.

```
stat_table <- getStatisticTable(ssc, ssc_params$params)
```
