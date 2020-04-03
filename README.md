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
Cardinal 2.4.0 or higher installed from Bioconductor 3.10 or higher.

## Usage

### spatialShrunkenCentroidsWrapper

Given multiple r, k, and s parameter sets, use this function to output one SpatialShrunkenCentroids2 object per set of parameters instead of one single SpatialShrunkenCentroids2 object containing every model.

```
ssc <- spatialShrunkenCentroids(data, r=c(1, 2), k=c(5, 10, 20), s=c(0, 3, 6, 9))

> class(ssc)
> SpatialShrunkenCentroids2
```
```
ssc <- spatialShrunkenCentroidsWrapper(data, r=c(1, 2), k=c(5, 10, 20), s=c(0, 3, 6, 9))

> class(ssc)
> List
```

### optimizeSSCParams

Use this function to assist in optimize the spatial neighborhood radius (r), initial number of segments (k), and sparsity (s) parameters used for running spatial shrunken centroids segmentation.

```
ssc <- spatialShrunkenCentroids(data, r=c(1, 2), k=c(5, 10, 20), s=c(0, 3, 6, 9))
optimizeSSCParams(ssc)
```

### updateMetadata

Use this function to add region of interest (ROI) and/or condition/treatment information based on x-y coordinates. An example Spot List (exported from Bruker FlexImaging v4.1) and ROI table can be found in the ```examples``` folder.

```
data <- updateMetadata(data, spotsFile="spots.txt", roiFile="rois.csv")
```

### spatialDGMMWrapper

Low quality features that can cause errors in spatialDGMM or cause the algorithm to get stuck when initializing Gaussian mixture models (GMMs). Use this function to apply spatial-DGMM segmentation on features one at a time and output a list of ```spatialDGMM``` objects. If a feature causes an error or is found to be unable to initalize GMMs, it will add a ```try-error``` object to the list at the index of said feature instead.

```
sdgmm <- spatialDGMMWrapper(data, r=1, k=4)
```

### segmentationTestWrapper

Use this function to following spatialDGMMWrapper to output a ```dataframe``` containing the results of ```segmentationTest()``` on each spatial-DGMM model. Only fixed conditions are currently supported. Since ```spatialDGMMWrapper()``` has each model contained in a separate ```spatialDGMM``` object, the AdjP values/FDR values in the resulting ```dataframe``` are unreliable and should be ignored.

```
sdgmm <- spatialDGMMWrapper(data, r=1, k=4)
resultsDf <- segmentationTestWrapper(sdgmmList, r=1, k=4, fixedCondition='treated', classControl='Ymax')
```

### findDiscriminantFeatures

Use this function to apply ```spatialDGMMWrapper()``` and ```segmentationTestWrapper()``` to an ```MSImagingExperiment``` object. Useful if only testing one fixed condition in the dataset.

```
resultsDf <- findDiscriminantFeatures(data, r=1, k=4, fixedCondition='treated', classControl='Ymax')
```
