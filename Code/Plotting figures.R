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
for (package in packages) {
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

## Create figures per each subset
path <- "Results"
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

for (f in folders) {
  
  pred_list <- list()
  
  response_dir <- file.path(f, "MaxEnt", "Response", "plots")
  
  png_files <- list.files(response_dir, pattern = "\\.png$", full.names = TRUE)
  png_files <- png_files[!grepl("ROC", basename(png_files), ignore.case = TRUE)]
  
  for (file_pred in png_files) {
    img <- image_read(file_pred)
    pred_list[[length(pred_list) + 1]] <- as.raster(image_convert(img, format = "rgba"))
  }
  
  grobs <- lapply(pred_list, rasterGrob, interpolate = TRUE)
  
  n <- length(grobs)
  nrow_grid <- 4
  ncol_grid <- 4
  
  layout_matrix <- matrix(
    1:(nrow_grid * ncol_grid),
    nrow = nrow_grid,
    ncol = ncol_grid,
    byrow = TRUE
  )
  
  tiff(file.path(f, "MaxEnt", "Response_curves.tiff"),
       width = 15 * 300,
       height = 14 * 300,
       res = 300,
       bg = "white")
  
  grid.arrange(
    grobs = grobs,
    layout_matrix = layout_matrix
  )
  
  dev.off()
}

## Create final figure
path <- "Results"
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

row_titles <- c(
  "General",
  "P1 Open-air funerary sites",
  "P1 Open-air settlemment sites",
  "P1 Cave sites",
  "P2 Open-air funerary sites",
  "P2 Open-air settlemment sites",
  "P2 Cave sites"
)

# Turn a PNG file name into a readable column title.
clean_var_name <- function(x) {
  x <- tools::file_path_sans_ext(basename(x))
  x <- gsub("_", " ", x)
  x
}

# Collect the response-curve PNGs (excluding ROC plots) for every folder,
# sorted so the same variable lands in the same column for every row.
png_by_folder <- lapply(folders, function(f) {
  response_dir <- file.path(f, "MaxEnt", "Response", "plots")
  png_files <- list.files(response_dir, pattern = "\\.png$", full.names = TRUE)
  png_files <- png_files[!grepl("ROC", basename(png_files), ignore.case = TRUE)]
  sort(png_files)
})

n_per_folder <- sapply(png_by_folder, length)
if (length(unique(n_per_folder)) != 1) {
  warning("Folders do not all contain the same number of response-curve PNGs; ",
          "columns may not line up correctly across rows.")
}
n_cols <- max(n_per_folder)

# Column titles taken from the folder with the most variables
ref_idx <- which.max(n_per_folder)
col_titles <- c(
  "Height",
  "TWI",
  "MSRM 3000",
  "Northing",
  "Easting",
  "Path frequency",
  "Path visibility (Top 10%)",  
  "Visibility Index",
  "Sky View Factor",
  "Agri suitability",
  "Cost from rivers and lakes",
  "Cost from coast",
  "Cost from variscite",
  "Cost from salt"
)

n_rows <- length(folders)

# Build a grob matrix: +1 row for column headers, +1 column for row headers
grob_matrix <- vector("list", (n_rows + 1) * (n_cols + 1))
dim(grob_matrix) <- c(n_rows + 1, n_cols + 1)

# Top-left corner: empty
grob_matrix[[1, 1]] <- nullGrob()

# Column headers (variable names)
for (j in seq_len(n_cols)) {
  grob_matrix[[1, j + 1]] <- textGrob(col_titles[j], gp = gpar(fontsize = 14, fontface = "bold"))
}

# Row headers (subset names) + the images themselves
for (i in seq_len(n_rows)) {
  grob_matrix[[i + 1, 1]] <- textGrob(row_titles[i], gp = gpar(fontsize = 14, fontface = "bold"), rot = 90)
  
  files_i <- png_by_folder[[i]]
  for (j in seq_len(n_cols)) {
    if (j <= length(files_i)) {
      img <- image_read(files_i[j])
      grob_matrix[[i + 1, j + 1]] <- rasterGrob(
        as.raster(image_convert(img, format = "rgba")),
        interpolate = TRUE
      )
    } else {
      grob_matrix[[i + 1, j + 1]] <- nullGrob()
    }
  }
}

# Flatten the matrix into the layout grid.arrange expects
grobs <- as.vector(t(grob_matrix))
layout_matrix <- matrix(seq_along(grobs), nrow = n_rows + 1, ncol = n_cols + 1, byrow = TRUE)

tiff(
  "Figure_all_response_curves_table.tiff",
  width  = (2 + 3 * n_cols) * 300,
  height = (1 + 3 * n_rows) * 300,
  res = 300,
  bg = "white"
)

grid.arrange(
  grobs = grobs,
  layout_matrix = layout_matrix,
  widths  = c(1.2, rep(3, n_cols)),
  heights = c(0.6, rep(3, n_rows))
)

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

# 4 Figure 7 and Supplementary 5.4 #############################################

path <- "Results"

# List all items in the directory that start with "p"
folders_starting_with_p <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with 'p' or 'P'
folders_starting_with_p <- folders_starting_with_p[!grepl("^z", basename(folders_starting_with_p), ignore.case = TRUE)]

# PDF filenames relative to folder
pdf_names <- c(
  "K_base.pdf",
  "Kres_pred_map.pdf",
  "KPPM_Ca.pdf",
  "KPPM_Ma.pdf",
  "KPPM_Th.pdf"
)

#Load the graphs
for (folder in folders_starting_with_p) {
  
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
  
  #Save the plot
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


# 5 Supplementary 4.4.1 ########################################################

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
