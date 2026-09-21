# README

This repository contains the code and data for:

**Investigating the influence of point placement strategies in the stability of modeled movement corridors: The case of the middle neolithic in the North-East of the Iberian Peninsula**

## Paper Authors

Biel Soriano Elias (a), Anna Bach Gómez (a) & Miquel Molist (a)

a Autonomous University of Barcelona, Prehistory Departament, SAPPO-GRAMPO


## Structure of the repository

The repository folder is structured as follows:

- **README.md**: This file (repository overview).  
- **Code/**: Contains all R scripts for the paper.  
  - **Main and functions**: Code (R scripts) to compute all the models in the paper, both functions and settings.
  - **Raster creation**:  Code (R scripts) to elaborate the experimental rasters used in the paper, not necessary to run as they are already provided.
- **Data/**: Raw data used in the paper, all in the CRS ETRS89 / UTM 31 N.  
  - **Rasters/**: All the nescessary rasters for the paper, including already computed experimental rasters and real world rasters
    - **Raster_a, Raster_b, Raster_c, Raster_dos_km & Raster_un_sis_kim**: All the experimental rasters created for the paper. The first three
    are intended for the first two steps of the experimental approach and the last two, altogether with Raster_a, for the last step.
    - **girona & girona_background**: The two needed rasters for the real world case of study of the paper. 
    Reprojected from the EU-DEM (https://ec.europa.eu/eurostat/web/gisco/geodata/digital-elevation-model/eu-dem, 
    last accessed on 9/1/2026).
  - **Vectors/**: All the nescessary vectors for the paper, including rivers, coast and sites.
    - **cat_ETRS89_31_N_modificat**: Shapefile (.shp and others) with the rivers used for the cost calculations in this paper. Note that some of them have
    been modified for the present analysis. Original data set was obtained from IGN (https://centrodedescargas.cnig.es/CentroDescargas/hidrografia, 
    last accessed on 9/1/2026)
    - **costa_girona**: Shapefile (.shp and others) with the coastline of the area under study. Data set obtained from ICGC 
    (https://www.icgc.cat/ca/Geoinformacio-i-mapes/Dades-i-productes/Geoinformacio-cartografica/Divisions-administratives, last accessed on 9/1/2026)
    - **punts_jaciments**: Shapefile (.shp and others) with sites used in the analysis.
    
## Computational Environment

All analyses and code development were conducted on:

- Windows 10 (64-bit): HP EliteBook 640 14inch G9 Notebook, 12th Gen Intel(R) Core (TM) i5, 16GB RAM  

### Information about the R Session

R version 4.5.2 (2025-10-31)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidyr_1.3.1            corrplot_0.95          stars_0.6-8            abind_1.4-8            spatstat_3.3-3         spatstat.linnet_3.2-6 
 [7] spatstat.model_3.3-6   rpart_4.1.24           spatstat.explore_3.4-3 nlme_3.1-168           spatstat.random_3.4-1  spatstat.geom_3.4-1   
[13] spatstat.univar_3.1-3  spatstat.data_3.1-6    dplyr_1.1.4            ggplot2_3.5.2          sf_1.0-21              terra_1.8-54          
[19] leastcostpath_2.0.12   raster_3.6-32          sp_2.2-0              

loaded via a namespace (and not attached):
 [1] gtable_0.3.6          xfun_0.52             spatstat.sparse_3.1-0 lattice_0.22-7        vctrs_0.6.5           tools_4.5.2          
 [7] spatstat.utils_3.1-4  generics_0.1.4        goftest_1.2-3         parallel_4.5.2        tibble_3.2.1          proxy_0.4-27         
[13] pkgconfig_2.0.3       Matrix_1.7-4          KernSmooth_2.23-26    RColorBrewer_1.1-3    lifecycle_1.0.4       compiler_4.5.2       
[19] farver_2.1.2          deldir_2.0-4          codetools_0.2-20      htmltools_0.5.8.1     class_7.3-23          yaml_2.3.10          
[25] pillar_1.10.2         classInt_0.4-11       iterators_1.0.14      foreach_1.5.2         tidyselect_1.2.1      digest_0.6.37        
[31] purrr_1.0.4           splines_4.5.2         polyclip_1.10-7       fastmap_1.2.0         grid_4.5.2            cli_3.6.5            
[37] magrittr_2.0.3        e1071_1.7-16          withr_3.0.2           tensor_1.5            scales_1.4.0          rmarkdown_2.29       
[43] evaluate_1.0.3        knitr_1.50            mgcv_1.9-3            rlang_1.1.6           Rcpp_1.0.14           glue_1.8.0           
[49] DBI_1.2.3             rstudioapi_0.17.1     R6_2.6.1              units_0.8-7 