# Cardinal Scripts

A small collection of scripts written to streamline data processing/analysis in Cardinal.

## Installation
#### Install from Github
devtools::install_github("gtluu/cardinalscripts")

#### Load from Local Folder
devtools::load_all("/path/to/cardinalscripts")

## Dependencies
Cardinal 2.4.0 or higher installed from Bioconductor 3.10 or higher.

## Usage

### optimizeSSCParams

Use this function to assist in optimize the spatial neighborhood radius (r), initial number of segments (k), and sparsity (s) parameters used for running spatial shrunken centroids segmentation.

### updateMetadata

Use this function to add region of interest (ROI) and/or condition/treatment information based on x-y coordinates. An example Spot List (exported from Bruker FlexImaging v4.1) and ROI table can be found in the ```examples``` folder.

### spatialShrunkenCentroidsWrapper

Given multiple r, k, and s parameter sets, use this function to output one SpatialShrunkenCentroids2 object per set of parameters instead of one single SpatialShrunkenCentroids2 object containing every model.
