##########################################################################
# This script is Variable creation for The settlement pattern            #
# of the Middle Neolithic of the North-eastern part of the Iberian       #   
# Peninsula: a homogeneous phenomenon?                                   #                                                   
# Biel Soriano-Elias, Anna Bach G?mez, Miquel Molist Montany?            #                                                   
#                                                                        #                                                   
# Author: Biel Soriano Elias                                             #
# Affiliation : Autonomous University of Barcelona                       #
# Creation date : 8/1/2026                                               #
# E-mail: biel.soriano@uab.cat                                           #
##########################################################################

# 0 Environment setup ##########################################################

# 0.1 Prepare environment ======================================================

# Set Working directory
setwd("XXXXX/XXXXX/XXXXX")

# Clean up workspace
rm(list=ls())

# 0.2 Install packages =========================================================

# Needed packages
packages <- c("leastcostpath","terra", "sf", "dplyr","tidyr", "doParallel", "viewscape",
              "ggplot2","spatstat","GA","stars","corrplot","maxnet","blockCV")

#Optional, Run this if the pacakges are not already installed
for (packages in packages) {
  install.packages(package, character.only = TRUE)
}

#Load the packages
for (package in packages) {
  library(package, character.only = TRUE)
}

# 0.3 Show session infos =======================================================

sessionInfo()

# 1 Creation fo variables ######################################################

# 1.0 DEM loading ==============================================================

#Loading the DEM 
rast_cat <- rast("Data/Rasters/11- Height.tiff")

# 2 Slope ######################################################################

#Compute Slope in degrees
pendent <- terrain(rast_cat, v="slope", unit="degrees", neighbors= 8)

# Save to file
writeRaster(pendent, "Data/Rasters/12- Slope in degrees.tiff", overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 3 Topographic Wetness Index (TWI) ############################################

# Compute flow direction
flow_dir <- terrain(rast_cat, v = "flowdir")

# Compute flow accumulation (upslope area)
flow_acc <- flowAccumulation(flow_dir)

# Flow accumulation (number of cells)
cell_area <- prod(res(rast_cat))
flow_acc_area <- flow_acc * cell_area

# Slope in radians
slope_rad <- terrain(rast_cat, v = "slope", unit = "radians")

# Compute TWI
twi <- log((flow_acc_area + 1) / tan(slope_rad))

#Erase Infinity results
twi[is.infinite(values(twi))] <- NA

# Save to file
writeRaster(twi, "Data/Rasters/13- TWI.tiff", overwrite = TRUE)


# 4 Multi-scale relief model ###################################################

## Function to msrm computation (Always forced to 1 fmin)
msrm <- function(r, fmin = 5, fmax = 100, x = 1.6, outdir = tempdir()) {
  # Compute resolution (m)
  rr <- mean(res(r))
  if (fmin <= rr) {
    message("fmin smaller than pixel size; setting fmin = raster resolution")
    fmin <- rr
  }
  
  # Calculate i and n values 
  i <- floor(((fmin - rr) / (2 * rr))^(1 / x))
  n <- ceiling(((fmax - rr) / (2 * rr))^(1 / x))
  
  # Ensure valid kernel indices
  if (is.nan(i) || i < 1) i <- 1
  if (is.nan(n) || n <= i) n <- i + 1
  
  message("Raster resolution: ", round(rr, 3), " m")
  message("i = ", i, ", n = ", n)
  
  reliefs <- c()
  prev <- NULL
  
  for (ndx in i:n) {
    rad <- round((ndx ^ x))
    if (rad < 1) rad <- 1
    if (rad > 50) {
      warning("Kernel radius > 50 pixels skipped to avoid memory issues.")
      next
    }
    
    k <- matrix(1, nrow = rad * 2 + 1, ncol = rad * 2 + 1)
    k <- k / sum(k)
    
    fpath <- file.path(outdir, paste0("LP_", ndx, ".tif"))
    f <- focal(r, w = k, fun = sum, na.rm = TRUE, filename = fpath, overwrite = TRUE)
    
    if (!is.null(prev)) {
      rpath <- file.path(outdir, paste0("RM_", ndx, ".tif"))
      rj <- prev - f
      writeRaster(rj, rpath, overwrite = TRUE)
      reliefs <- c(reliefs, rpath)
    }
    prev <- f
  }
  
  # Combine and compute mean of relief models
  msrm_files <- rast(reliefs)
  msrm_raw <- mean(msrm_files)
  msrm <- round(msrm_raw * 1000) / 1000
  
  return(msrm)
}  

## Compute MSRM
msrm_result <- msrm(rast_cat, fmin = 300, fmax = 3000, x = 1)

# Save to file
writeRaster(msrm_result, "Data/Rasters/14- MSRM.tiff", overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 5 Northing and Easting  ######################################################

#Compute aspect
aspect_raster <- terrain(rast_cat, v = "aspect")

#Convert aspect from degrees to radians
aspect_radians <- aspect_raster * pi / 180

##Calculate the continuous East and North components
#Easting (X) = sin(aspect in radians)
easting_raster <- sin(aspect_radians)

#Northing (Y) = cos(aspect in radians)
northing_raster <- cos(aspect_radians)

#Save the results to files
writeRaster(easting_raster, "Data/Rasters/16- Easting.tiff")
writeRaster(northing_raster, "Data/Rasters/15- Northing.tif")

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 6 Path Visibility  ###########################################################

### Caution, computationally expensive ###

#Load the necessary data
y <- st_read("Data/Vectors/rast_cat_100_points.shp")
z <- rast("Data/Rasters/21- Path frequency.tiff")
elevacio <- rast_cat

#Set up storing folders
dir.create("Data/Rasters/Visibility")
dir.create("Data/Rasters/Visibility/RASTERS_PARTS_VC")
dir.create("Data/Rasters/Visibility/RASTERS_VC")

# Obtaining raster extent
ext_r <- ext(elevacio)
xmin <- ext_r[1]
xmax <- ext_r[2]
ymin <- ext_r[3]
ymax <- ext_r[4]

# Defining block size in cells
block_cells <- 200  # Number of cells per block in each dimension
res_x <- res(elevacio)[1]  # Resolution in x direction
res_y <- res(elevacio)[2]  # Resolution in y direction

# Compute block size in spatial units
block_size_x <- block_cells * res_x
block_size_y <- block_cells * res_y

# List to store the blocks
blocks <- list()
index <- 1

# Divide the raster into blocks
for (i in seq(xmin, xmax, by = block_size_x)) {
  for (j in seq(ymin, ymax, by = block_size_y)) {
    
    # Ensure last blocks fit within raster by adjusting size dynamically
    x_end <- min(i + block_size_x, xmax)
    y_end <- min(j + block_size_y, ymax)
    
    # Defining block extension
    ext_block <- ext(i, x_end, j, y_end)
    
    # Crop the raster block
    block_raster <- crop(elevacio, ext_block)
    
    # Write raster block to file
    is_empty <- all(is.na(values(block_raster))) # avoid saving empty rasters
    
    # Creation of the blocks
    if (is_empty) {
      
      print("Skipping raster: It is entirely NA")
      
    } else {
      
      blocks[[index]] <- block_raster
      writeRaster(blocks[[index]], filename = file.path("Data/Rasters/Visibility/RASTERS_PARTS_VC", paste0("block_", index, ".tif")), overwrite = TRUE)
      
      print(paste("Block", index, "created with extension", ext_block))
      index <- index + 1
    }
  }
}

rm(blocks)

gc()

print("Step 1 completed")

ncores <- detectCores() - 1 #identification of number of pc cores

# List all TIFF files in the RASTERS_RAW folder
llista_rast_parts <- list.files("Data/Rasters/Visibility/RASTERS_PARTS_VC", pattern = "\\.tif[f]?$", full.names = TRUE)

numeros <- as.numeric(gsub("Data/Rasters/Visibility/RASTERS_PARTS_VC/block_|\\.tif", "", llista_rast_parts))

# Order raster by number 
llista_rast_parts_sorted <- llista_rast_parts[order(numeros)]

# Read all rasters
for(r in seq_along(llista_rast_parts_sorted)) {
  
  r <- r
  
  # Loading the raster's block
  elevacio_r <- rast(llista_rast_parts_sorted[r])
  
  ext_rast <- ext(elevacio_r)
  res_rast <- res(elevacio_r)
  raster_buit <- rast(ext = ext_rast, res = res_rast, vals = NA, crs = crs(elevacio))
  
  polygons_ras <- as.polygons(elevacio_r, dissolve = TRUE)
  names(polygons_ras)[1] <- "Layer_1" 
  polygons_ras_fil <- polygons_ras[polygons_ras$Layer_1 >= 0, ]
  polygons_ras_fil_sf <- st_as_sf(aggregate(polygons_ras_fil))
  buff_3000 <- st_buffer(polygons_ras_fil_sf, 3000) # for paths creation
  elevacio_r_3000 <- crop(elevacio, buff_3000)
  
  rm(elevacio_r)
  rm(polygons_ras)
  rm(polygons_ras_fil)
  rm(buff_3000)
  
  punts_r <- st_intersection(y, polygons_ras_fil_sf)
  
  rm(polygons_ras_fil_sf)
  
  viewshed <- compute_viewshed(elevacio_r_3000, punts_r, r = 3000, parallel = TRUE)
  rast_r_3000 <- crop(z, elevacio_r_3000)
  rast_r_3000 <- resample(rast_r_3000, elevacio_r_3000)
  
  rm(elevacio_r_3000)
  
  gc()
  
  #Loop to calculate the viewshed of each point
  for (i2 in 1:nrow(punts_r)) {
    
    #Superposing viewshed raster and path raster
    viewshed_punt <- viewshed[[as.character(i2)]] #extracting viewshed
    rast_view <- visualize_viewshed(viewshed_punt, outputtype = "raster") #converting it to raster
    crop_rast <- crop(rast_r_3000, rast_view)
    extracted_raster <- mask(crop_rast, rast_view, maskvalues = 1, inverse = TRUE) #masking path raster with viewshed raster
    
    # Extract values from path raster
    extracted_values <- terra::values(extracted_raster, ID = FALSE)
    
    # Calculate the 90th quantile
    punts_r[i2,"VC"] <- quantile(extracted_values, 0.9, na.rm = TRUE)
    
  }
  
  rm(viewshed)
  
  raster_buit <- rasterize(punts_r, raster_buit, field = "VC")
  
  #saving of the raster's file  
  nom_arxiu <- paste("block",as.character(r),"total_viewshed_paths.tif") #creation of file name
  writeRaster(raster_buit, filename = file.path("Data/Rasters/Visibility/RASTERS_VC",nom_arxiu)) #saving of the raster of viewshed
  
  rm(punts_r)
  rm(raster_buit)
  
  gc()
  
}


# List all TIFF files in the RASTERS_VC folder
llista_rast <- list.files(file.path("Data/Rasters/Visibility/RASTERS_VC"), pattern = "\\.tif[f]?$", full.names = TRUE)

# Read all rasters
rast_list <- lapply(llista_rast, terra::rast)

# Read all rasters
rast_list_norm <- lapply(rast_list, terra::rast)
rast_list_norm <- sprc(rast_list_norm)

# Build the mosaic
mosaic_rast <- mosaic(rast_list_norm, fun = mean)

writeRaster(mosaic_rast, filename = "Data/Rasters/22- Path Visibility (Top 10%).tif") #saving of the raster of viewshed

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 7 Visibility Index  ##########################################################

#Load the necessary data
punts <- st_read("Data/Vectors/rast_cat_100_points.shp")
elevacio <- rast_cat

#Set up storing folders
dir.create("Data/Rasters/Visibility/RASTERS")
dir.create("Data/Rasters/Visibility/RASTER_PARTS")


# Obtaining raster extent
ext_r <- ext(elevacio)
xmin <- ext_r[1]
xmax <- ext_r[2]
ymin <- ext_r[3]
ymax <- ext_r[4]

# Defining block size in cells
block_cells <- 200  # Number of cells per block in each dimension
res_x <- res(elevacio)[1]  # Resolution in x direction
res_y <- res(elevacio)[2]  # Resolution in y direction

# Compute block size in spatial units
block_size_x <- block_cells * res_x
block_size_y <- block_cells * res_y

# List to store the blocks
blocks <- list()
index <- 1

# Divide the raster into blocks
for (i in seq(xmin, xmax, by = block_size_x)) {
  for (j in seq(ymin, ymax, by = block_size_y)) {
    
    # Ensure last blocks fit within raster by adjusting size dynamically
    x_end <- min(i + block_size_x, xmax)
    y_end <- min(j + block_size_y, ymax)
    
    # Defining block extension
    ext_block <- ext(i, x_end, j, y_end)
    
    # Crop the raster block
    block_raster <- crop(elevacio, ext_block)
    
    # Write raster block to file
    is_empty <- all(is.na(values(block_raster))) # avoid saving empty rasters
    
    # Creation of the blocks
    if (is_empty) {
      
      print("Skipping raster: It is entirely NA")
      
    } else {
      
      blocks[[index]] <- block_raster
      writeRaster(blocks[[index]], filename = file.path("Data/Rasters/Visibility/RASTER_PARTS", paste0("block_", index, ".tif")), overwrite = TRUE)
      
      print(paste("Block", index, "created with extension", ext_block))
      index <- index + 1
    }
  }
}

rm(blocks)

print("Step 1 completed")

ncores <- detectCores() - 1 #identification of number of pc cores

# List all TIFF files in the RASTERS_RAW folder
llista_rast_parts <- list.files("Data/Rasters/Visibility/RASTER_PARTS", pattern = "\\.tif[f]?$", full.names = TRUE)

numeros <- as.numeric(gsub("Data/Rasters/Visibility/RASTER_PARTS/block_|\\.tif", "", llista_rast_parts))

# Oder rasters by number
llista_rast_parts_sorted <- llista_rast_parts[order(numeros)]


# Read all rasters
for(r in seq_along(llista_rast_parts_sorted)) {
  
  r <- r
  
  # Loading the raster's block
  elevacio_r <- rast(llista_rast_parts_sorted[r])
  
  ext_rast <- ext(elevacio_r)
  res_rast <- res(elevacio_r)
  raster_buit <- rast(ext = ext_rast, res = res_rast, vals = NA, crs = crs(elevacio))
  
  polygons_ras <- as.polygons(elevacio_r, dissolve = TRUE)
  names(polygons_ras)[1] <- "Layer_1" 
  polygons_ras_fil <- polygons_ras[polygons_ras$Layer_1 >= 0, ]
  polygons_ras_fil_sf <- st_as_sf(aggregate(polygons_ras_fil))
  buff_3000 <- st_buffer(polygons_ras_fil_sf, 3000) # for paths creation
  elevacio_r_3000 <- crop(elevacio, buff_3000)
  
  rm(elevacio_r)
  rm(polygons_ras)
  rm(polygons_ras_fil)
  rm(buff_3000)
  
  punts_r <- st_intersection(punts, polygons_ras_fil_sf)
  
  rm(polygons_ras_fil_sf)
  
  viewshed <- compute_viewshed(elevacio_r_3000, punts_r, r = 3000, parallel = TRUE, workers = ncores)
  
  rm(elevacio_r_3000)
  
  gc()
  
  #Loop to calculate the viewshed of each point
  for (i2 in 1:nrow(punts_r)) {
    
    punts_r[i2, "cells"] <- sum(viewshed[[as.character(i2)]]@visible)
    
    if (i2 == 1) {rast_max <- ncell(viewshed[["1"]]@visible)}
    
  }
  
  rm(viewshed)
  
  raster_buit <- rasterize(punts_r, raster_buit, field = "cells")
  
  rast_min <- 0
  
  raster_buit <- ((raster_buit - rast_min)/(rast_max - rast_min)) 
  
  #saving of the raster's file  
  nom_arxiu <- paste("block",as.character(r),"total_viewshed.tif") #creation of file name
  writeRaster(raster_buit, filename = file.path("Data/Rasters/Visibility/RASTERS",nom_arxiu)) #saving of the raster of viewshed
  
  rm(punts_r)
  rm(raster_buit)
  
  gc()
  
}


# List all TIFF files in the RASTERS_RAW folder
llista_rast <- list.files(file.path("Data/Rasters/Visibility/RASTERS"), pattern = "\\.tif[f]?$", full.names = TRUE)

# Read all rasters
rast_list <- lapply(llista_rast, terra::rast)
rast_list <- sprc(rast_list)

# Build the mosaic (may have to be done in smaller parts for proper working)
mosaic_rast <- mosaic(rast_list, fun = mean)

#Save the results
writeRaster(mosaic_rast, filename = "Data/Rasters/23- Visibility Index.tif") 

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 8 Visual Prominence Index ####################################################

#Loading total viewshed raster
total_view <- rast("Data/Rasters/23- Visibility Index.tiff")

#Settign up MSRM function
msrm <- function(r, fmin = 5, fmax = 100, x = 1.6, outdir = tempdir()) {
  # --- Compute resolution (m) ---
  rr <- mean(res(r))
  if (fmin <= rr) {
    message("fmin smaller than pixel size; setting fmin = raster resolution")
    fmin <- rr
  }
  
  # --- Calculate i and n values ---
  i <- floor(((fmin - rr) / (2 * rr))^(1 / x))
  n <- ceiling(((fmax - rr) / (2 * rr))^(1 / x))
  
  # Ensure valid kernel indices
  if (is.nan(i) || i < 1) i <- 1
  if (is.nan(n) || n <= i) n <- i + 1
  
  message("Raster resolution: ", round(rr, 3), " m")
  message("i = ", i, ", n = ", n)
  
  reliefs <- c()
  prev <- NULL
  
  for (ndx in i:n) {
    rad <- round((ndx ^ x))
    if (rad < 1) rad <- 1
    if (rad > 50) {
      warning("Kernel radius > 50 pixels skipped to avoid memory issues.")
      next
    }
    
    k <- matrix(1, nrow = rad * 2 + 1, ncol = rad * 2 + 1)
    k <- k / sum(k)
    
    fpath <- file.path(outdir, paste0("LP_", ndx, ".tif"))
    f <- focal(r, w = k, fun = sum, na.rm = TRUE, filename = fpath, overwrite = TRUE)
    
    if (!is.null(prev)) {
      rpath <- file.path(outdir, paste0("RM_", ndx, ".tif"))
      rj <- prev - f
      writeRaster(rj, rpath, overwrite = TRUE)
      reliefs <- c(reliefs, rpath)
    }
    prev <- f
  }
  
  # Combine and compute mean of relief models
  msrm_files <- rast(reliefs)
  msrm_raw <- mean(msrm_files)
  msrm <- round(msrm_raw * 1000) / 1000
  
  return(msrm)
}

## Visual prominence (difference of a cell in regard neighbors in total viewshed values) with MSRM (Orengo i Petrie, 2018) 
visual_prominance_msrm <-  msrm(total_view, fmin = 300, fmax = 3000, x = 1) 

# Normalize the raster to the range [0, 1]
min_value_2 <- minmax(visual_prominance_msrm)[1]
max_value_2 <- minmax(visual_prominance_msrm)[2]

visual_normalized <- (visual_prominance_msrm - min_value_2) / (max_value_2 - min_value_2)

#Save the results
writeRaster(visual_normalized, filename = "Data/Rasters/25- Visual Prominance index 3000 (MSRM).tif") 

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 9 Agricultural suitability  ##################################################

# 9.1 Extracting reference data of rast_cat ====================================

#Obtaining most of the data from terra package
ref_raster <- rast_cat
bb <- ext(ref_raster)
in_crs <- crs(ref_raster) #CRS of the area of interest


dir.create("Data/Rasters/Soil")
dir.create("Data/Rasters/Soil/RAW")
dir.create("Data/Rasters/Soil/Class")

# 9.2 PH Data ==================================================================

#Output setting
resolution = "250m" 

#Variables selection
voi = "phh2o" # variable of interest
depth = "15-30cm"
quantile = "mean"
voi_layer = paste(voi,depth,quantile, sep="_") # layer of interest
rst_crs = "ESRI:54052"

rstFile = paste0("/vsicurl/https://files.isric.org/soilgrids/latest/data/", voi, '/', voi_layer,'.vrt')

rst = rast(rstFile)
crs(rst) = rst_crs # assign crs as the ESRI code rather than proj string

bb_proj = project(bb,from=crs(in_crs),to=crs(rst_crs)) #project bb to same crs as raster layer

window(rst) = NULL # remove any existing window
window(rst) = bb_proj # get just roi

#Saving the result
writeRaster(rst, paste0("Data/Rasters/Soil/RAW/","PH_250.tiff"), overwrite = TRUE)


# 9.3 Type of soils data =======================================================

voi = "wrb" 
soil_class = "MostProbable" 

voi_layer = soil_class

rst_crs = "EPSG:4326" # EPSG code for the WRB layers

rstFile = paste0("/vsicurl/https://files.isric.org/soilgrids/latest/data/", voi, '/', voi_layer,'.vrt')

rst = rast(rstFile)
crs(rst) = rst_crs # assign crs as the ESRI code rather than proj string

bb_proj = project(bb,from=crs(in_crs),to=crs(rst_crs)) #project bb to same crs as raster layer

window(rst) = NULL # remove any existing window
window(rst) = bb_proj # get just roi

#Saving the result
writeRaster(rst, paste0("Data/Rasters/Soil/RAW/","Soil_type_250.tiff"), overwrite = TRUE)

# 9.4 Soil depth (200) data ====================================================

rst_crs = "EPSG:4326" # EPSG code for the WRB layers

rstFile = "/vsicurl/https://files.isric.org/soilgrids/former/2017-03-10/data/BDRICM_M_250m_ll.tif"

rst = rast(rstFile)
crs(rst) = rst_crs # assign crs as the ESRI code rather than proj string

bb_proj = project(bb,from=crs(in_crs),to=crs(rst_crs)) #project bb to same crs as raster layer

window(rst) = NULL # remove any existing window
window(rst) = bb_proj # get just roi

#saving the results
writeRaster(rst, paste0("Data/Rasters/Soil/RAW/","BD_200_250.tiff"), overwrite = TRUE)

# 9.5 Soil depth (total) data ==================================================

rst_crs = "EPSG:4326" # EPSG code for the WRB layers

rstFile = "/vsicurl/https://files.isric.org/soilgrids/former/2017-03-10/data/BDTICM_M_250m_ll.tif"

rst = rast(rstFile)
crs(rst) = rst_crs # assign crs as the ESRI code rather than proj string

bb_proj = project(bb,from=crs(in_crs),to=crs(rst_crs)) #project bb to same crs as raster layer

window(rst) = NULL # remove any existing window
window(rst) = bb_proj # get just roi

#Saving the results
writeRaster(rst, paste0("Data/Rasters/Soil/RAW/","BD_full_250.tiff"), overwrite = TRUE)


# 9.6 Reclasifiyng rasters and creating the index raster =======================

### Reprojecting and extending each raster to match reference
raster_files <- list.files("Data/Rasters/Soil/RAW", pattern = "\\.tif[f]?$", full.names = TRUE)
raster_list <- lapply(raster_files, rast) 

rast_names <- c("BD_200_100",
                "BD_full_100",
                "PH_100",
                "Soil_type_100")

raster_list <- lapply(seq_along(raster_list), function(i) {
  
  r <- raster_list[[i]]
  
  # Reproject to reference CRS
  if (crs(r) != crs(ref_raster)) {
    if (i %in% c(4)) {  
      r <- project(r, ref_raster, method = "near")  # categorical
    } else { 
      r <- project(r, ref_raster, method = "bilinear")  # continuous
    }
  }
  
  # Extend to reference extent
  if (!all(ext(r) == ext(ref_raster))) {
    r <- extend(r, ext(ref_raster))
  }
  
  # Match resolution
  if (!all(res(r) == res(ref_raster))) {
    if (i %in% c(4)) {
      r <- resample(r, ref_raster, method = "near")
    } else {
      r <- resample(r, ref_raster, method = "bilinear")
    }
  }
  
  # Apply mask using the reference raster
  r <- mask(r, ref_raster)
  
  # Save the raster
  writeRaster(r, filename = file.path("Data/Rasters/Soil", paste0(rast_names[i], ".tiff")), overwrite = TRUE)
  
  return(r)
  
})

### Reclassifying rasters for index construction
#Loading reprojected rasters
raster_files <- list.files("Data/Rasters/Soil", pattern = "\\.tif[f]?$", full.names = TRUE)
raster_list <- lapply(raster_files, rast) 

## PH (look up table from Salvador Baiges, 2024 - p. 162 - Blat)
rast_ph <- raster_list[[3]]

#Filling the gaps
filled <- focal(rast_ph, w = matrix(1, 51, 51), fun = mean, na.policy = "only", na.rm = TRUE)
filled_final <- mask(filled, ref_raster)
writeRaster(filled_final, "Data/Rasters/Soil/PH_filled_100.tiff")
raster_list[[3]] <- rast("Data/Rasters/Soil/PH_filled_100.tiff")

#Reclassify the raster
classify_ph <- function(x) {
  out <- rep(NA, length(x))
  out[x >= 65 & x <= 75] <- 5
  out[(x >= 60 & x < 65) | (x > 75 & x <= 82)] <- 4
  out[(x >= 56 & x < 60) | (x > 82 & x <= 83)] <- 3
  out[(x >= 52 & x < 56) | (x > 83 & x <= 85)] <- 2
  out[(x >= 42 & x < 52) | (x > 85 & x <= 86)] <- 1
  out[x < 42 | x > 86] <- 0
  return(out)
}

classified_ph <- app(rast_ph, classify_ph)
writeRaster(classified_ph, "Data/Rasters/Soil/Class/PH_class_100.tiff")

## BD (No full as look up table does not support it - Salvador Baiges, 2024 - p. 163 - Ordi)
rast_bd_200 <- raster_list[[1]]

classify_cm <- function(x) {
  out <- rep(NA, length(x))
  out[x > 150] <- 5
  out[x >= 100 & x <= 150] <- 4
  out[x >= 50 & x < 100] <- 3
  out[x >= 20 & x < 50] <- 1
  out[x < 20] <- 0
  return(out)
}

classified_bd_200 <- app(rast_bd_200, classify_cm)
writeRaster(classified_bd_200, "Data/Rasters/Soil/Class/BD_200_class_100.tiff")

## Soil type (look up table from Salvador Baiges, 2024 - p. 164 - Intensiu/Planes al·luvials)
rast_st <- raster_list[[4]]

classify_soil_ref <- function(x) {
  out <- rep(0, length(x))  # default class 0
  out[x >= 25 & x <= 35] <- 5
  out[x >= 76 & x <= 84] <- 5
  out[x >= 48 & x <= 52] <- 5
  out[x >= 66 & x <= 67] <- 4
  out[x >= 21 & x <= 24] <- 3
  out[x >= 68 & x <= 72] <- 2
  out[x >= 97 & x <= 98] <- 1
  return(out)
}

classified_st <- app(rast_st, classify_soil_ref)

#mask values inside the polygon
classified_st <- mask(classified_st, ref_raster)
writeRaster(classified_st, "Data/Rasters/Soil/Class/Soil_class_100.tiff")

## Slope (look up table from  Salvador Baiges, 2024 - p. 161 - Intensiu/Planes al·luvials)
rast_slope <- rast("Data/Rasters/12- Slope in degrees.tiff")

classify_percentage <- function(x) {
  out <- rep(NA, length(x))
  out[x >= 0 & x <= 5] <- 5
  out[x > 5 & x <= 8] <- 4
  out[x > 8 & x <= 16] <- 3
  out[x > 16 & x <= 30] <- 1
  out[x > 30] <- 0
  return(out)
}

classified_slope <- app(rast_slope, classify_percentage)
writeRaster(classified_slope, "Data/Rasters/Soil/Class/Slope_class_100.tiff")

### Combining the rasters into one single layer
#Load the classified rasters in a list
raster_files <- list.files("Data/Rasters/Soil/Class", pattern = "\\.tif[f]?$", full.names = TRUE)
raster_list <- lapply(raster_files, rast) 
raster_stack <- rast(raster_list)

#Sum
sum_raster <- app(raster_stack, sum, na.rm = TRUE)

#Normalize between 0 and 1
min_val <- global(sum_raster, "min", na.rm = TRUE)[1,1]
max_val <- global(sum_raster, "max", na.rm = TRUE)[1,1]

normalized_raster <- (sum_raster - min_val) / (max_val - min_val)

#Save the resutls
writeRaster(normalized_raster, "Data/Rasters/31- Agri suitability.tif", overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat")))
gc()

# 10 Cost from rivers, coast, salt and variscite ###############################

# 10.1 Aggregate DEM ===========================================================

#Aggregate DEM to facilitate cost computation
r_lr <- aggregate(rast_cat, fact= 2, fun = mean)
cost_raster <- create_slope_cs(r_lr, cost_function =  "tobler", neighbours = 16)

# 10.2 Cost from rivers and lakes ==============================================

#Load the unified rivers and lakes
rivers <- st_read("Data/Vectors/Rius_&_llacs.shp")

#Sample points every 30 m along each line
samples <- st_line_sample(rivers, density = 1 / 200, type = "regular")  

#Initiate an empty object to store cleaned points
all_points <- st_sfc(crs = st_crs(rivers))

#Clean points
for (i in seq_along(samples)) {
  if (length(samples[[i]]) > 0) {
    # Wrap as sfc before casting — preserves CRS
    mp <- st_sfc(samples[[i]], crs = st_crs(rivers))
    pts <- st_cast(mp, "POINT")
    all_points <- c(all_points, pts)
  }
}

rivers_points <- st_sf(geometry = all_points)  #Convert to sf object

# Initialize a raster to store the minimum accumulated cost values
min_cc <- rasterise(cost_raster)
values(min_cc) <- Inf

#Compute cost for all the point
for (i in 1:nrow(rivers_points)) {
  # Compute accumulated cost from each origin
  coords <- st_coordinates(rivers_points[i,])
  
  # Check if point is inside raster extent
  if (!all(coords[,1] >= xmin(min_cc) & coords[,1] <= xmax(min_cc) &
           coords[,2] >= ymin(min_cc) & coords[,2] <= ymax(min_cc))) {
    next  # skip this iteration
  }
  
  cc <- create_accum_cost(x = cost_raster, origins = rivers_points[i,], FUN = mean, rescale = FALSE)
  
  # Update min_cc with the minimum value between the existing and the new cc
  min_cc <- min(min_cc, cc, na.rm = TRUE)
  
  rm(cc)
  gc()
}

# Replace Inf with NA
min_cc[values(min_cc) == Inf] <- NA

# Save to file
writeRaster(min_cc, "Data/Rasters/32- Cost from rivers and lakes", overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat","cost_raster","r_lr")))
gc()

# 10.3 Cost from coast =========================================================

#Load the coast 
coast <- st_read("Data/Vectors/Rast_cat_coast_line_mod.shp")

#Sample points every 30 m along each line
samples <- st_line_sample(coast, density = 1 / 200, type = "regular")  

#Initiate an empty object to store cleaned points
all_points <- st_sfc(crs = st_crs(coast))

#Clean points
for (i in seq_along(samples)) {
  if (length(samples[[i]]) > 0) {
    # Wrap as sfc before casting — preserves CRS
    mp <- st_sfc(samples[[i]], crs = st_crs(coast))
    pts <- st_cast(mp, "POINT")
    all_points <- c(all_points, pts)
  }
}

coast_points <- st_sf(geometry = all_points)  #Convert to sf object

# Initialize a raster to store the minimum accumulated cost values
min_cc <- rasterise(cost_raster)
values(min_cc) <- Inf

#Compute cost for all the point
for (i in 1:nrow(coast_points)) {
  # Compute accumulated cost from each origin
  coords <- st_coordinates(coast_points[i,])
  
  # Check if point is inside raster extent
  if (!all(coords[,1] >= xmin(min_cc) & coords[,1] <= xmax(min_cc) &
           coords[,2] >= ymin(min_cc) & coords[,2] <= ymax(min_cc))) {
    next  # skip this iteration
  }
  
  cc <- create_accum_cost(x = cost_raster, origins = coast_points[i,], FUN = mean, rescale = FALSE)
  
  # Update min_cc with the minimum value between the existing and the new cc
  min_cc <- min(min_cc, cc, na.rm = TRUE)
  
  rm(cc)
  gc()
}

# Replace Inf with NA
min_cc[values(min_cc) == Inf] <- NA

# Save to file
writeRaster(min_cc, "Data/Rasters/33- Cost from coast", overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat","cost_raster","r_lr")))
gc()

# 10.4 Cost from salt ==========================================================

#Load points
cambrils <- c(X=367020.9,Y=4665977.1) #Salí de cambrils
cardona <- c(X=390228.5,Y=4640144.6) #Muntanya de sal de Cardona
salt <- rbind(cambrils,cardona)

## Compute cost
#Compute raster
min_cc <- rasterise(cost_raster)
values(min_cc) <- NA

#Compute cost of each point
for (i in 1:nrow(salt)) {
  
  cc <- create_accum_cost(x = cost_raster, origins = salt[i,], FUN = mean, rescale = FALSE)
  
  # Update min_cc with the minimum value between the existing and the new cc
  min_cc <- min(min_cc, cc, na.rm = TRUE)
  
  rm(cc)
  rm(coords)
  gc()
}

#Replace Inf with NA
min_cc[values(min_cc) == Inf] <- NA
min_cc <- resample(min_cc, rast_cat, method = "bilinear")

#Adapt to final res and ext
ref_raster <- rast_cat
ref_extent <- ext(ref_raster)
ref_crs <- crs(ref_raster)
ref_res <- res(ref_raster)

# Reproject to reference CRS
if (crs(min_cc) != crs(ref_raster)) {
  min_cc <- project(min_cc, ref_raster, method = "bilinear")  # continuous
}

# Extend to reference extent
if (!all(ext(min_cc) == ext(ref_raster))) {
  min_cc <- extend(min_cc, ext(ref_raster))
}

# Match resolution
if (!all(res(min_cc) == res(ref_raster))) {
  min_cc <- resample(min_cc, ref_raster, method = "bilinear")
}


# Apply mask using the reference raster
min_cc <- mask(min_cc, rast_cat)

# Save the results
writeRaster(min_cc, file.path("Data/Rasters", "35- Cost from salt.tif"), overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat","cost_raster","r_lr")))
gc()

# 10.5 Cost from variscite =====================================================

#Load point
general <- st_read("Conjunts/general/general.shp")
mina <- general[general$FID_ == 1.269 ,]

#Compute cost from point
cc <- create_accum_cost(x = cost_raster, origins = mina, FUN = mean, rescale = FALSE)
cc <- resample(cc, rast_cat, method = "bilinear")

#Adapt to final res and ext
ref_raster <- rast_cat
ref_extent <- ext(ref_raster)
ref_crs <- crs(ref_raster)
ref_res <- res(ref_raster)

# Reproject to reference CRS
if (crs(cc) != crs(ref_raster)) {
  cc <- project(cc, ref_raster, method = "bilinear")  # continuous
}

# Extend to reference extent
if (!all(ext(cc) == ext(ref_raster))) {
  cc <- extend(cc, ext(ref_raster))
}

# Match resolution
if (!all(res(cc) == res(ref_raster))) {
  cc <- resample(cc, ref_raster, method = "bilinear")
}

# Apply mask using the reference raster
cc <- mask(cc, ref_raster)

# Save the results
writeRaster(cc, file.path("Data/Rasters", "34- Cost from variscita.tif"), overwrite = TRUE)

rm(list=setdiff(ls(), c("rast_cat","cost_raster","r_lr")))
gc()