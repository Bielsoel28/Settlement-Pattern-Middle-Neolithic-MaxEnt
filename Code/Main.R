##########################################################################
# This script is Main for The settlement pattern of the Middle Neolithic #
# of the North-eastern part of the Iberian Peninsula:                    #   
# a homogeneous phenomenon?                                              #                                                   
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
packages <- c("terra", "sf", "dplyr","tidyr",
              "ggplot2","corrplot","SDMtune","blockCV",
              "nortest","rJava","doParallel")

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

# 1 Creating subsets of sites ##################################################

##Loading sites points
loc <-  st_read("Data/Vectors/sites.shp")
loc <- loc[!st_is_empty(loc), ]

#Function to keep only unique sites
keep_first_unique <- function(df, column_name) {
  df[!duplicated(df[[column_name]]), ]
}

## Creating subsets of sites
#General 
loc_general <- keep_first_unique(loc,"Name")

#NMI
loc_p1 <- loc[loc$Period == "NMI", ]

##Cave
loc_p1_cov <- loc_p1[grepl("\\(Cave\\)", loc_p1$Type, ignore.case = TRUE), ]

##Open-air Funerary and habitat sites
loc_p1_all <- loc_p1[!grepl("\\(Cave\\)", loc_p1$Type, ignore.case = TRUE), ]
loc_p1_all_f <- loc_p1_all[grepl("Sepulture", loc_p1_all$Type, ignore.case = TRUE), ]
loc_p1_all_h <- loc_p1_all[grepl("Settlement", loc_p1_all$Type, ignore.case = TRUE), ]

rm(loc_p1_all)
rm(loc_p1)

#NMP
loc_p2 <- loc[loc$Period == "NMP", ]

##Cave
loc_p2_cov <- loc_p2[grepl("\\(Cave\\)", loc_p2$Type, ignore.case = TRUE), ]

##Open-air Funerary and habitat sites
loc_p2_all <- loc_p2[!grepl("\\(Cave\\)", loc_p2$Type, ignore.case = TRUE), ]
loc_p2_all_f <- loc_p2_all[grepl("Sepulture", loc_p2_all$Type, ignore.case = TRUE), ]
loc_p2_all_h <- loc_p2_all[grepl("Settlement", loc_p2_all$Type, ignore.case = TRUE), ]

rm(loc_p2_all)
rm(loc_p2)

#Creating a list of point sets for loop functionality
list_punts <- list(loc_general, 
                   loc_p1_all_f,loc_p1_all_h, 
                   loc_p1_cov, 
                   loc_p2_all_f,loc_p2_all_h, 
                   loc_p2_cov)
#Freeing memory
rm(loc_general, 
   loc_p1_all_f,loc_p1_all_h, 
   loc_p1_cov, 
   loc_p2_all_f,loc_p2_all_h, 
   loc_p2_cov)

gc()

#Creating a vector with names for readability of the results
names(list_punts) <- c("general",
                       "p1_all_f", "p1_all_h", 
                       "p1_cov",
                       "p2_all_f", "p2_all_h",
                       "p2_cov")

## Creating a folders to store them
#Creating Main folder
dir.create("Results")

#Creating a folder for each of the subsets
invisible(lapply(names(list_punts), function(name) {
  dir_path <- file.path("Results", name)
  
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  } else {
    message("Skipping: ", dir_path, " - directory already exists.")
  }
}))

#Saving each of them
sapply(seq_along(list_punts), function(i) {
  
  obj_name <- names(list_punts)[i]
  
  # Build file path
  file_path <- file.path(
    "Results",
    obj_name,
    paste0(obj_name, ".shp")
  )
  
  # If file already exists, skip writing
  if (file.exists(file_path)) {
    message(paste("Skipping:", file_path, "- file already exists."))
    return(NULL)
  }
  
  # Otherwise write the spatial file
  st_write(list_punts[[i]], file_path, quiet = TRUE)
})

# 2 Normality and correlation of variables #####################################

# 2.1 Setting up data ==========================================================

## creating folder to store results
dir.create("Results/Variables_cor")

## Loading rasters
raster_files <- list.files("Data/Rasters", pattern = "\\.tif[f]?$", full.names = TRUE)
raster_list <- lapply(raster_files, rast) 

#Stack rasters
predictors <- rast(raster_list)

#Assign names
names(predictors) <- tools::file_path_sans_ext(basename(raster_files))

### Extracting values of the variables
#Loading mesh of background points
sample <- st_read("Data/Vectors/fishnet_rast_cat_100_red_points_v3.shp")
#Extracting values
sample_pred <- terra::extract(predictors, sample, ID = FALSE) 

#Erasing the NA points from sample for background aggregation correct functionality
sample <- cbind(sample_pred, sample)
sample <- na.omit(sample)
sample <- st_as_sf(sample)

#Saving final values
sample_pred <- na.omit(sample_pred)

# 2.2 Anderson-Darling Normality test ==========================================

## Loop for AD test
normality_results <- lapply(sample_pred, function(column) {
  if (is.numeric(column) && length(column) >= 3) {
    ad.test(column)$p.value
  } else {
    NA
  }
})

# Convert results to a data frame
normality_results_df <- data.frame(
  Column = names(normality_results),
  P_Value = unlist(normality_results)
)

# Add a column to interpret the results
normality_results_df$Normal <- ifelse(normality_results_df$P_Value > 0.05, "Yes", "No")

# Save the results
write.csv(normality_results_df, paste0("Results/Variables_cor","Normality.csv"))

# 2.3 Pearson correlation test test ============================================

## Pearson test
sample_pred_cor <- cor(sample_pred, method = "pearson")

#Ploting and saving the result
tiff(file.path("Results/Variables_cor","Correlation_plot_pearson.tiff"), width = 9*300, height = 6*300, res = 300) # Width and height in pixels
corrplot.mixed(sample_pred_cor,lower.col = "black", number.cex = 0.35, tl.pos="lt", tl.cex= 0.4, cl.cex = 0.4)
dev.off()

## Filtering out correlated variables

predictors_final <- predictors[[-c(2,11)]]
sample_values <- sample_pred[,-c(2,11)]

rm(list=setdiff(ls(),c("list_punts","predictors_final","sample","sample_values")))

gc()

# 3 Selection of background points #############################################
### Folllowing background aggregation method (Xu et al., 2024)

# 3.1 Setting up data for parallel computing ===================================

## Decay factor for geographic distance
k <- 1  

## Loading Presence points
presences <- list_punts[[1]]  

# Extract coordinates for background and presence points
bg_coords <- st_coordinates(sample)
pres_coords <- st_coordinates(presences)

# Extract environmental values for the  presence points
pres_values <- terra::extract(predictors_final, presences, ID = FALSE)

##Prepare objects for loop functionality
n <- nrow(presences)
aggregation_scores <- numeric(nrow(sample))  # Vector to store aggregation scores

# 3.2 Loop for parallel calculation aggregation scores =========================

## Set up parallel loop 
ncores <- 10
sim_per_core <- ceiling(nrow(sample) / ncores)  # split iterations across cores

# Create cluster
clus <- makeCluster(ncores)

# Load required packages on each worker
clusterEvalQ(clus, {
  library(terra)
  library(sf)
})

# Export necessary objects to workers
clusterExport(clus, varlist = c("sample", "bg_coords", "sample_values", "pres_coords", "pres_values", "k", "n", "aggregation_scores"))

## Worker function
worker_aggregation <- function(i) {
  # Extract coordinates for background and presence points
  bg_coords_i <- bg_coords[i, ]
  bg_values_i <- sample_values[i, ]
  
  # Calculate geographic distance (d_j) between background point and all presence points
  distances <- sqrt((bg_coords_i[1] - pres_coords[, 1])^2 + (bg_coords_i[2] - pres_coords[, 2])^2 + 1^2)
  
  # Calculate weight (w_j) based on geographic distance (d_j)
  weights <- 1 / (distances^k)
  
  # Calculate the standard deviation (sigma_v) for the environmental variable values at presence points
  sigma_v <- apply(pres_values, 2, sd)
  
  # Calculate the bandwidth (h) using Equation (2)
  h <- sigma_v * (4 / (3 * n))^0.2
  
  # Initialize environmental similarity (S_v_i)
  Sv_i <- 0
  
  # Loop over all presence points (j) to compute the weighted sum of environmental similarities
  for (j in 1:n) {
    diff_env <- (bg_values_i - pres_values[j, ])
    similarity_term <- exp(- (diff_env^2) / (2 * h^2))
    Sv_i <- Sv_i + weights[j] * (1 / sqrt(2 * pi)) * similarity_term
  }
  
  # Normalize Sv_i by n * h (as per the original formula)
  Sv_i <- Sv_i / (n * h)
  
  # Return the aggregation score for the current background point
  return(mean(as.numeric(Sv_i)))
}

# Run the function in parallel for each row of the sample
aggregation_scores <- unlist(parLapply(clus, 1:nrow(sample), worker_aggregation))

# Stop the cluster after computations are done
stopCluster(clus)

# Normalize the aggregation scores
normalized_scores <- (aggregation_scores - min(aggregation_scores)) / (max(aggregation_scores) - min(aggregation_scores))

## Sample background points based on the normalized aggregation scores
set.seed(234)  # Set seed for reproducibility
sampled_background_ids <- sample(sample$FID_fishne, size = 10000, prob = normalized_scores, replace = FALSE)

# Select the sampled background points from the original dataset
sample_bg_ag_sampled <- sample[sample$FID_fishne %in% sampled_background_ids, ]
sample_bg_ag <- sample
sample_bg_ag$Normalized_scor <- normalized_scores #save all values

# Plot the results
tiff(file.path("Results/Variables_cor","selected_bg_points.tiff"), width = 9*300, height = 6*300, res = 300) # Width and height in pixels
plot(st_geometry(sample), col = "gray", pch = 16, main = "Presence, Background, and Sampled Background Points", 
     xlab = "Longitude", ylab = "Latitude", cex = 0.6, axes = TRUE)

plot(st_geometry(sample_bg_ag_sampled), col = "blue", pch = 20, add = TRUE, cex = 0.8)

plot(st_geometry(presences), col = "red", pch = 20, add = TRUE, cex = 0.8)

legend("bottomright", 
       legend = c("Background Points", "Presence Points", "Sampled Background Points"), 
       col = c("gray", "red", "blue"), 
       pch = c(16, 20, 20),  # Different shapes for different point types
       pt.cex = c(0.6, 0.8, 0.8),  # Adjust point size for better visibility
       bty = "n",  # Remove the border around the legend
       cex = 0.9)  # Adjust font size

dev.off()

# Save the selected background points 
st_write(sample_bg_ag, file.path("Results/Variables_cor","all_bg_points.shp"))
st_write(sample_bg_ag_sampled, file.path("Results/Variables_cor","selected_bg_points.shp"))

rm(list=setdiff(ls(),c("list_punts","predictors_final","sample","sample_bg_ag_sampled")))

gc()


# 4 MaxEnt modelling ###########################################################

# 4.1 Cross-folds validation blocks creation ===================================

tiff(file.path("Results/Variables_cor","Spatial_blocks_red_mod_final.tiff"), width = 12*300, height = 8*300, res = 300) # Width and height in pixels
sac <- cv_spatial_autocor(predictors_final) 
dev.off()

#Load DEM for plotting
rast_cat <- rast("Data/Rasters/11- Height.tiff")

# 4.2 MaxEnt computing =========================================================

## Loop for sites subsets
for (i in names(list_punts)){ 
  
  #Folder creation for storing results
  out_dir <- file.path("Results",i,"MaxEnt")
  dir.create(out_dir) #General folder
  
  #Extracting loc values
  loc <- list_punts[[i]]
  jac_values <- terra::extract(predictors_final, loc, ID = FALSE)
  
  # Preparing the data
  jac_total <- cbind(jac_values, loc$geometry)
  sample_values <- terra::extract(predictors_final, sample_bg_ag_sampled, ID = FALSE) 
  sample_total <- cbind(sample_values, sample_bg_ag_sampled$geometry)
  sample_total <- na.omit(sample_total)
  
  jac_total$Id <- 1
  sample_total$Id <- 0
  
  all_data <- data.frame(rbind(jac_total, sample_total))

  ##Loop to ensure valid kfolds
  set.seed(123) #set seed
  valid_split <- FALSE
  
  #Main loop
  set.seed(123) #set seed
  valid_split <- FALSE
  iteration_count <- 0
  max_iterations <- 10
  num_folds <- 5  # start with 5 folds
  
  while (!valid_split) {
    
    kfolds <- cv_spatial(
      x = all_data,
      column = "Id",
      r = predictors_final,
      size = sac$range, 
      k = num_folds,
      selection = "random",
      iteration = 100,
      hexagon = FALSE,
      plot = FALSE
    )
    
    # extract number of test points per fold
    counts <- kfolds[["records"]][[paste0("test_", 1)]]
    
    # check the condition
    if (all(counts >= 3)) {
      valid_split <- TRUE
    } else {
      iteration_count <- iteration_count + 1
      
      # if 10 attempts fail and still invalid, switch to 3 folds
      if (iteration_count >= max_iterations && num_folds == 5) {
        num_folds <- 3
        iteration_count <- 0  # reset counter for 3-fold attempt
      }
    }
  }
  
  ## Saving and plotting Blokcs results
  # Spatial blocks without points
  gg <- cv_plot(kfolds, r = rast_cat)  # returns a ggplot object
  ggsave(filename = file.path(out_dir, "Spatial_blocks_mde.tiff"),
         plot = gg, width = 8, height = 6, dpi = 300, units = "in")
  
  # Spatial blocks with sampling points
  gg2 <- cv_plot(kfolds, x = all_data, r = rast_cat)
  ggsave(filename = file.path(out_dir, "Spatial_blocks_points_mde.tiff"),
         plot = gg2, width = 8, height = 6, dpi = 300, units = "in")
  
  ## Creating a SWD object with data
  data_swd <- new("SWD",
                  species = i,
                  coords  = as.data.frame(st_coordinates(all_data$geometry)),
                  data    = all_data[, 1:ncol(sample_values)],
                  pa      = as.numeric(all_data$Id))  
    
    
    #GA optimization
    first_model <- train("Maxent", 
                         data = data_swd,
                         folds = kfolds,
                         progress = TRUE)
    
    h <- list(fc = c("l", "lq", "lh", "lqp", "lqph", "lqpht"),
              reg = seq(0.2, 5, 0.2),
              iter = c(500,1000,2000))
    
    
    genetic <- optimizeModel(first_model, 
                             hypers = h, 
                             metric = "auc", 
                             pop = 20, 
                             gen = 5,
                             keep_best = 0.4,
                             keep_random = 0.2,
                             mutation_chance = 0.4,
                             interactive = FALSE,
                             progress = TRUE)
    
    best_id <- which.max(genetic@results$test_AUC)
    best_fc   <- genetic@results$fc[best_id]
    best_reg  <- genetic@results$reg[best_id]
    best_iter <- genetic@results$iter[best_id]
    
    MaxEnt_model <- train(
      "Maxent",
      data  = data_swd,
      fc    = best_fc,
      reg   = best_reg,
      iter  = best_iter,
      folds = kfolds,
      progress = TRUE
    )
    
    # Return results 
    results <- list(
      MaxEnt_model = MaxEnt_model,
      train_data = data_swd,
      kfolds = kfolds
    )
  
  save(results, file = file.path(out_dir,paste0("MaxEnt_models.Rdata")))
  
}

rm(list=setdiff(ls(),c("list_punts","predictors_final")))

gc()

# 4.3 Model evaluation =========================================================

auc_df <- data.frame(
  AUC = numeric(length(names(list_punts))),
  row.names = names(list_punts)
)

# Load the functions
source("./Code/Auxiliar/plotROC_kfold.R")
source("./Code/Auxiliar/modelReportCV.R")

##Loop to run over all models
for (i in names(list_punts)){ 
  
##Set folder where to retrieve models
out_dir <- file.path("Results",i,"MaxEnt")

#Load the models
load(file.path(out_dir,paste0("MaxEnt_models.Rdata")))

##AUC 
auc_df[i,] <- SDMtune::auc(results$MaxEnt_model, test = TRUE)

## Response curves
modelReportCV(model = results$MaxEnt_model,
              folder = file.path(out_dir,"Response"),
              type = "cloglog",
              response_curves = TRUE,
              only_presence = FALSE,
              permut = 10,
              verbose = TRUE)

##Variable importance
var <- varImp(results$MaxEnt_model, permut = 10, progress = FALSE)

#Plot the data
var_plot <- plotVarImp(var, color = "#abd9e9")
ggsave(file.path(out_dir,paste0("Variable importance.tiff")), var_plot, width = 6, height = 8)

#save the data
write.csv(var, file = file.path(out_dir,"Permutation_Importance.csv"))

##Jackknife test
jk <- doJk(results$MaxEnt_model,
           metric = "auc",
           with_only = TRUE,
           progress = FALSE)

#Plot the data
jK <- plotJk(jk, type = c("train"))

ggsave(file.path(out_dir,paste0("Jackknife_test.tiff")), jK, width = 6, height = 8)

#Save the data
write.csv(jk, file = file.path(out_dir,"Jackknife_test.csv"))

##Prediction map creation
names(predictors_final) <- names(results$MaxEnt_model@data@data)

#Pred map from every folds
pred_list <- lapply(
  results$MaxEnt_model@models,
  function(m) {
    terra::predict(
      predictors_final,
      m,
      clamp = FALSE,
      type = "cloglog",
      na.rm = TRUE
    )
  }
)

#Mean of all folds
pred_stack <- rast(pred_list)
mean_pred_map <- app(pred_stack, mean, na.rm = TRUE)

writeRaster(mean_pred_map, file.path(out_dir, "pred_map.tiff"), overwrite = TRUE)

}

write.csv(auc_df, file = "Results/AUC_ME.csv", row.names = TRUE)

