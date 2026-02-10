# README

This repository contains the code and data for:

**The settlement pattern of the Middle Neolithic of the North-eastern part of the Iberian Peninsula: a homogeneous phenomenon?**

## Paper Authors

Biel Soriano Elias (a), Anna Bach Gómez (a) & Miquel Molist (a)

a Autonomous University of Barcelona, Prehistory Departament, SAPPO-GRAMPO

## Structure of the repository

The repository folder is structured as follows:

- **README.md**: This file (repository overview).  
- **Code/**: Contains all R scripts for the paper.  
  - **Main**: Code (R script) to perform MaxEnt modelling and all associated analysis.
  - **Plotting figures**: Code (R script) to plot Fig. 2, Fig.5, Fig. 6, Fig. 7 and Supplementay 6 figures.
  - **Point Process models**: Code (R script) to perform and fitt all PPMs of the analysis
  - **Variables creation**: Code (R script) to create all the varaibles used in the nalysis, not nescessary to run as they all are provided.
  - **Auxiliars/**: Contains auxiliar scripts for functions used in R scripts for the paper
    - **modelReportCV**: Code (R script) to create reports of MaxEnt models
    - **modelReportCV.html**: Template for MaxEnt model report
    - **plotROC_kfold**: Code (R script) to plot AUC graphics for MaxEnt models
- **Data/**: Data used in the paper, all in the CRS ETRS89 / UTM 31 N.  
  - **Rasters/**: All the nescessary rasters for the paper
    - **11- Height**: Resampled version (100 x 100 meters) of the original 25 x 25 meter Digital Elevation Model (DEM) raster supplied by GLO-30 Copernicus 
    (https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/copernicus#Elevation, last accessed on 15/12/2025 at 15:55)
    - **12- Slope in degrees**: Slope on degrees raster
    - **13- TWI**: Topographic Wetness Index raster
    - **14- MSRM 3000**: MSRM raster
    - **15- Northing**: Northing raster
    - **16- Easting**: Easting raster
    - **21- Path frequency**: Path Frequency raster
    - **22- Path visibility (Top 10%)**: Path Visibility raster
    - **23- Visibility Index**: Total viewshed Index raster
    - **24- Sky View Factor**: Sky View Factor raster 
    - **25- Visual Prominance index 3000 (MSRM)**: Visual Prominance Index raster
    - **31- Agri suitability**: Agricultural suitability index raster
    - **32- Cost from rivers and lakes**: Cost from rivers and lakes raster
    - **33- Cost from coast**: Cost from coast raster
    - **34- Cost from variscita**:  Cost from coast variscite mines
    - **35- Cost from salt**: Cost from slat outcrops
  - **Vectors/**: All the necessary vectors for the paper 
    - **Buff_27450_mod**: Shapefile (.shp and others) with the area under study limits
    - **fishnet_rast_cat_100_red_points_v3**: Shapefile (.shp and others) with mesh of 1km background points
    - **rast_cat_100_points**: Shapefile (compressed in .zip) with points from the vectorization of "11- Height" for visibility variables construction
    - **Rast_cat_coast_line_mod**: Shapefile (.shp and others) with modified coastline for "33- Cost from coast", original data set obtained from EEA
    (https://www.eea.europa.eu/en/datahub/datahubitem-view/af40333f-9e94-4926-a4f0-0a787f1d2b8f, last accessed on 9/1/2025)
    - **Rius_&_llacs**: Shapefile (.shp and others) with modified Rivers and Lakes, original data sets obtained from IGN  
    (https://centrodedescargas.cnig.es/CentroDescargas/hidrografia, last accessed on 9/1/2026) and ICGC 
    (https://www.icgc.cat/ca/Geoinformacio-i-mapes/Mapes/Mapa-de-cobertes-del-sol-de-Catalunya, last accessed on 9/1/2026)
    - **sites**: Shapefile (.shp and others) with the sites used in the analysis
  - **Sites references**: All the references for the database of sites (.pdf)
- **Results/**: Folder to store the results of the analysis, its creation is also scripted in "Main" 
  - **Variables_cor/**: Folder to store results of correlation analysis of variables, its creation is also scripted in "Main" 
    - **selected_bg_points**: Shapefile (.shp and others) with selected background points for the MaxEnt models
  
## Computational Environment

All analyses and code development were conducted on:

- Windows 10 (64-bit): HP EliteBook 640 14inch G9 Notebook, 12th Gen Intel(R) Core (TM) i5, 16GB RAM  

### Information about the R Session

R version 4.5.2 (2025-10-31 ucrt)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rJava_1.0-11  nortest_1.0-4 blockCV_3.1-5 SDMtune_1.3.2 corrplot_0.95 ggplot2_3.5.2 tidyr_1.3.1   dplyr_1.1.4  
 [9] sf_1.0-21     terra_1.8-54 

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       compiler_4.5.2     tidyselect_1.2.1   Rcpp_1.0.14        scales_1.4.0       lattice_0.22-7    
 [7] R6_2.6.1           generics_0.1.4     classInt_0.4-11    dismo_1.3-16       tibble_3.2.1       units_0.8-7       
[13] DBI_1.2.3          pillar_1.10.2      RColorBrewer_1.1-3 rlang_1.1.6        sp_2.2-0           cli_3.6.5         
[19] withr_3.0.2        magrittr_2.0.3     class_7.3-23       grid_4.5.2         rstudioapi_0.17.1  lifecycle_1.0.4   
[25] vctrs_0.6.5        KernSmooth_2.23-26 proxy_0.4-27       glue_1.8.0         raster_3.6-32      farver_2.1.2      
[31] codetools_0.2-20   e1071_1.7-16       purrr_1.0.4        tools_4.5.2        pkgconfig_2.0.3   