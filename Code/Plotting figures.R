##########################################################################
# This script is Plotting figures for The settlement pattern             #
# of the Middle Neolithic of the North-eastern part of the Iberian       #   
# Peninsula: a homogeneous phenomenon?                                   #                                                   
# Biel Soriano-Elias, Anna Bach G?mez, Miquel Molist Montany?            #                                                   
#                                                                        #                                                   
# Author: Biel Soriano Elias                                             #
# Affiliation : Autonomous University of Barcelona                       #
# Creation date : 8/1/2025                                               #
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

# Load extra packages
library(magick)
library(grid)
library(gridExtra)

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

# Load the grpahs
for (f in seq_along(folders)) {
  file_pred <- file.path(folders[[f]], "MaxEnt", "Combined_Response_Curves_All_Variables.tiff")
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
  grid.text(clean_titles[i],
            vp = viewport(layout.pos.row=row_col[1], layout.pos.col=row_col[2]),
            y = unit(1, "npc") - unit(2, "lines"),
            gp=gpar(fontsize=14,fontface="bold"))
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
  file_pred <- file.path(folder_path, "pred_map_total.tiff")
  
  if (!file.exists(file_pred)) {
    message("Skipping ", folder, ": pred_map_total.tiff not found.")
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
  file.path("Figure response.tiff"),
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

#Load extra packages
library(grid)
library(gridExtra)
library(pdftools)
library(magick)

# Select folders and files 
path <- "Results"

# List all items in the directory 
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with "v"
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

# PDF filenames relative to folder
pdf_names <- c(
  "K_base.pdf",
  "Kres_pred_map.pdf",
  "KPPM_Ca.pdf",
  "KPPM_Ma.pdf",
  "KPPM_Th.pdf"
)

## Loop to go over all folders
for (folder in folders) {
  
  folder_path <- file.path(folder, "MaxEnt")
  
  pdf_grobs <- list()
  for (pdf_name in pdf_names) {
    pdf_file <- file.path(folder_path, pdf_name)
    if (!file.exists(pdf_file)) next
    
    img <- image_read_pdf(pdf_file, density = 300)
    img <- img[1]  # first page
    
    r <- as.raster(image_convert(img, format = "rgba"))
    pdf_grobs[[length(pdf_grobs)+1]] <- rasterGrob(r, interpolate = TRUE)
  }
  
  if (length(pdf_grobs) == 0) next
  
  num_pdfs <- length(pdf_grobs)
  
  # Determine layout: 2 columns per row
  n_col <- 2
  n_row <- ceiling(num_pdfs / n_col)
  
  layout_matrix <- matrix(seq_len(n_row * n_col), nrow = n_row, ncol = n_col, byrow = TRUE)
  
  # Replace empty slots with NA if num_pdfs < n_row * n_col
  layout_matrix[layout_matrix > num_pdfs] <- NA
  
  # Output TIFF
  out_file <- file.path(folder_path, "Combined_PPM.tiff")
  
  tiff(
    out_file,
    width = 8*300,
    height = 12*300,
    res = 300,
    bg = "white"
  )
  
  grid.arrange(
    grobs = pdf_grobs,
    layout_matrix = layout_matrix
  )
  
  dev.off()
  message("Saved combined PDFs for folder: ", folder)
}