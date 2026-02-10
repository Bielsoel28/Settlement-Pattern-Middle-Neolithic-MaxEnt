##########################################################################
# This script is Plotting figures for The settlement pattern             #
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
packages <- c("terra", "sf", "dplyr","tidyr","ggplot2","magick","grid","gridExtra","pdftools","gtable")

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

# 1 Plotting figure 2 ##########################################################

#Load the rasters
raster_files <- list.files("Data/Rasters", pattern = "\\.tif[f]?$", full.names = TRUE)
raster_list <- lapply(raster_files, rast) 
predictors <- rast(raster_list)
names(predictors) <- tools::file_path_sans_ext(basename(raster_files))

#Plot them
tiff(
  file.path("Figure_variables.tiff"),
  width  = 12 * 300,
  height = 9 * 300,
  res    = 300
)

# number of layers
n <- nlyr(predictors)

# rows needed for max 3 per row
par(mfrow = c(ceiling(n / 3), 3), mar = c(2, 2, 2, 4))

plot(predictors)

dev.off()

# 2 Plotting figure 5 ##########################################################

## Collect each subset plots
path <- "Results"
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

for (f in folders) {
  
  pred_list <- list()
  
  response_dir <- file.path(f, "MaxEnt", "Response","plots")
  
  # List all PNG files
  png_files <- list.files(response_dir, pattern = "\\.png$", full.names = TRUE)
  
  # Remove files containing "roc" (case-insensitive)
  png_files <- png_files[!grepl("ROC", basename(png_files), ignore.case = TRUE)]
  
  # Read each PNG file
  for (file_pred in png_files) {
    img <- image_read(file_pred)
    pred_list[[length(pred_list) + 1]] <- as.raster(image_convert(img, format = "rgba"))
  }
  
  # Convert to grobs
  grobs <- lapply(pred_list, rasterGrob, interpolate=TRUE)
  
  # Layout: 4 x 4
  n <- length(grobs)
  nrow_grid <- 4
  ncol_grid <- 4
  
  # If fewer images than grid spaces, fill with NULLs
  layout_matrix <- matrix(
    1:(nrow_grid * ncol_grid),
    nrow = nrow_grid,
    ncol = ncol_grid,
    byrow = TRUE
  )
  
  # Save final composed figure
  tiff(file.path(f, "MaxEnt","Response_curves.tiff"),
       width = 15*300,
       height = 14*300,
       res = 300,
       bg = "white")
  
  grid.arrange(
    grobs = grobs,
    layout_matrix = layout_matrix
  )
  
  
  dev.off()
  
}

## Create the final figure
# Select folders and files 
path <- "Results"

# List all items in the directory 
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with "v"
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

# Titles
clean_titles <- c(
  "General",
  "P1 Open-air funerary sites",
  "P1 Open-air settlemment sites",
  "P1 Cave sites",
  "P2 Open-air funerary sites",
  "P2 Open-air settlemment sites",
  "P2 Cave sites"
)

## Load images
#Create an empty object for loading the graphs
pred_list <- list()

# Load the graphs
for (f in seq_along(folders)) {
  file_pred <- file.path(folders[[f]], "MaxEnt", "Response_curves.tiff")
  if (!file.exists(file_pred)) next
  img <- image_read(file_pred)
  if (!is.null(img)) {
    pred_list[[length(pred_list)+1]] <- as.raster(image_convert(img, format = "rgba"))
  }
}

# Convert to grobs
grobs <- lapply(pred_list, rasterGrob, interpolate=TRUE)

# Layout: 1 big top + 2-column bottom
n <- length(grobs)
bottom_indices <- if (n > 1) 2:n else c()
layout_matrix <- rbind(
  c(1,1),
  if(length(bottom_indices) > 0) matrix(bottom_indices, ncol = 2, byrow = FALSE)
)

# Plot the graphs
tiff(
  "Figure_all_response_curves.tiff",
  width = 18*300,
  height = 30*300,
  res = 300,
  bg = "white"
)

grid.arrange(
  grobs = grobs,
  layout_matrix = layout_matrix,
  top = textGrob(clean_titles[1], gp=gpar(fontsize=16,fontface="bold"))
)

# Add titles for small plots manually
pushViewport(viewport(layout=grid.layout(nrow(layout_matrix), ncol(layout_matrix))))
for (i in seq_along(grobs)) {
  row_col <- which(layout_matrix == i, arr.ind = TRUE)
  if (i == 1) next
  grid.text(
    clean_titles[i],
    vp = viewport(layout.pos.row = row_col[1],
                  layout.pos.col = row_col[2]),
    y = unit(0.95, "npc") + unit(0.5, "lines"),
    just = "bottom",
    gp = gpar(fontsize = 14, fontface = "bold")
  )
}

dev.off()

# 3 Plotting figure 6 ##########################################################

# Select folders and files 
path <- "Results"

# List all items in the directory 
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with "v"
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

##Load all rasters in a single raster 
pred_list <- list()

for (f in seq_along(folders)) {
  
  folder <- folders[[f]]
  folder_path <- file.path(folder, "MaxEnt")
  file_pred <- file.path(folder_path, "pred_map.tiff")
  
  if (!file.exists(file_pred)) {
    message("Skipping ", folder, ": pred_map.tiff not found.")
    next
  }
  
  r <- rast(file_pred)
  
  pred_list[[folder]] <- r
}

# Combine into a single SpatRaster
pred_raster <- rast(pred_list)

## Plot them 
n <- nlyr(pred_raster)

zlim <- range(values(pred_raster), na.rm = TRUE)

cols <- hcl.colors(100, "Red-Green", rev = TRUE)

clean_titles <- c("General","P1 Open-air funerary sites", "P1 Open-air settlemment sites", "P1 Cave sites",
                  "P2 Open-air funerary sites", "P2 Open-air settlemment sites", "P2 Cave sites")

tiff(
  file.path("Figure pred maps.tiff"),
  width  = 9 * 300,
  height = 15 * 300,
  res    = 300
)
n <- nlyr(pred_raster)

# Number of bottom plots
n_bottom <- n - 1

# Number of rows needed for 2 columns
n_rows <- ceiling(n_bottom / 2)

# Build layout matrix (by column)
bottom_mat <- matrix(
  seq_len(n_bottom) + 1,
  nrow = n_rows,
  ncol = 2,
  byrow = FALSE
)

# Add the big top map
layout_mat <- rbind(
  c(1, 1),
  bottom_mat
)

layout(
  layout_mat,
  heights = c(2.2, rep(1.7, nrow(layout_mat) - 1))
)

par(mar = c(1, 1, 1, 1))

plot(
  pred_raster[[1]],
  col = cols,
  zlim = zlim,
  legend = FALSE,
  main = clean_titles[1]
)

# Legend inside top map
plot(
  pred_raster[[1]],
  col = cols,
  zlim = zlim,
  legend.only = TRUE,
  smallplot = c(0.88, 0.93, 0.15, 0.85),
  axis.args = list(
    at = pretty(zlim, 5),
    labels = round(pretty(zlim, 5), 2)
  ),
  legend.args = list(
    text = "Prediction value",
    side = 4,
    line = 2,
    cex = 0.9
  )
)

# Remaining maps
for (i in 2:n) {
  plot(
    pred_raster[[i]],
    col = cols,
    zlim = zlim,
    legend = FALSE,
    main = clean_titles[i]
  )
}


dev.off()

# 4 Plotting figure 7 and supplementary 6 ######################################

# Select folders and files 
path <- "Results"

# List all items in the directory 
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with "v"
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

# Titles
clean_titles <- c(
  "General",
  "P1 Open-air funerary sites",
  "P1 Open-air settlemment sites",
  "P1 Cave sites",
  "P2 Open-air funerary sites",
  "P2 Open-air settlemment sites",
  "P2 Cave sites"
)

## Load images
#Create an empty object for loading the graphs
pred_list <- list()

# Load the graphs
for (f in seq_along(folders)) {
  file_pred <- file.path(folders[[f]], "MaxEnt", "Spatial_blocks_points_mde.tiff")
  if (!file.exists(file_pred)) next
  img <- image_read(file_pred)
  if (!is.null(img)) {
    pred_list[[length(pred_list)+1]] <- as.raster(image_convert(img, format = "rgba"))
  }
}

# Convert to grobs
grobs <- lapply(pred_list, rasterGrob, interpolate=TRUE)

# Layout: 1 big top + 2-column bottom
n <- length(grobs)
bottom_indices <- if (n > 1) 2:n else c()
layout_matrix <- rbind(
  c(1,1),
  if(length(bottom_indices) > 0) matrix(bottom_indices, ncol = 2, byrow = FALSE)
)

# Plot the graphs
tiff(
  "Figure_Spatial_blocks.tiff",
  width = 18*300,
  height = 30*300,
  res = 300,
  bg = "white"
)

grid.arrange(
  grobs = grobs,
  layout_matrix = layout_matrix,
  top = textGrob(clean_titles[1], gp=gpar(fontsize=16,fontface="bold"))
)

# Add titles for small plots manually
pushViewport(viewport(layout=grid.layout(nrow(layout_matrix), ncol(layout_matrix))))
for (i in seq_along(grobs)) {
  row_col <- which(layout_matrix == i, arr.ind = TRUE)
  if (i == 1) next
  grid.text(clean_titles[i],
            vp = viewport(layout.pos.row=row_col[1], layout.pos.col=row_col[2]),
            y = unit(1, "npc") - unit(2, "lines"),
            gp=gpar(fontsize=14,fontface="bold"))
}

dev.off()