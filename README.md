# Nonlinear Trends and Regional Climate Signals in Lake Ice Cover Duration

## Abstract

Lake ice phenology is a sensitive indicator of climate change, however, trends in ice cover duration vary significantly by region. This study examines long-term changes in lake ice cover duration using singular spectrum analysis (SSA) to identify temporal trends and periodic patterns.  Satellite-derived brightness temperature was incorporated as a proxy for surface air temperature to investigate climatic drivers of phenological change. Results reveal a significant decline in ice cover duration in Sweden and Finland, while comparable trends were not detected in other examined regions. A significant correlation between brightness temperature and ice cover duration was observed in Sweden and Finland, but not elsewhere. Spatial correlation analysis further revealed stronger negative relationships between brightness temperature and ice cover duration in the Nordic region, whereas correlations in the Northeastern United States were weaker and more variable. These findings suggest that regions characterized by lower mean temperatures and moderate seasonal variability may exhibit more pronounced changes in lake ice phenology. By integrating time series decomposition with multivariate analysis, this study provides a statistical framework for identifying regional differences in lake ice responses to changing climatic conditions. 

## Contents

This repository contains the R scripts and relevant data files downloaded from the [NSIDC]{https://doi.org/10.7265/N5W66HP8} and [NOAA]{https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ncdc:C01578/html} for recreating all the results from this project. The contents are as follows:

- `data\`
  - `liag_ice.csv`: data from the NSIDC of ice cover durations of lakes and rivers in the Northern Hemisphere.
  - `liag_info.csv`: data from the NSIDC of physical characteristics of those lakes and rivers.
  - `uat4_tb_v04r00_avrg_chtlt_s197812_e202312_c20240104.nc`: global data from the NOAA of monthly average 2.5 degree gridded brightness temperatures of the Lower Troposphere.
  - `corrs20.csv`: data of the correlation between yearly ice cover duration and the yearly average brightness temperature.
- `Rscripts/`
  - `icessacode.R`: general R script containing code for SSA and generation of `corrs20.csv`. 
  - `spatial_interpolation.R`: R script for inverse-distance weighting. 

## Download

Downloading the data files are unnecessary, as the R scripts pull the data from this repository themselves. However, due to the nature of the brightness temperature data file type, the data set must be downloaded to be read into R. This is done in the script `icessacode.R`, though it requires the `icessacode.R` script to be saved locally. The script then checks the source file location for the existence of the brightness temperature data, and if it does not find it downloads the file from this repository. So, *make sure you save* `icessacode.R` *before running the code to ensure the data is pulled correctly*. 


