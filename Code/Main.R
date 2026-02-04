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
packages <- c("leastcostpath","terra", "sf", "dplyr","tidyr", "doParallel", "viewscape",
              "ggplot2","spatstat","GA","stars","corrplot","maxnet","blockCV","nortest")

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

# 4.2 Functions for MaxEnt evaluation =====================================

#Function to compute Jackknife test
jackknife_test <- function(results_list, 
                           response_col = "presence", 
                           type = "cloglog") {
  # Initialize list for fold-wise jackknife
  jack_list <- lapply(results_list, function(res) {
    if(!res$skipped) {
      model <- res$MaxEnt_model
      data  <- res$train_data
      
      vars <- colnames(data)
      vars <- vars[vars != response_col]
      
      y_true <- data[[response_col]]
      
      # Compute AUC with each variable alone and without each variable
      fold_res <- lapply(vars, function(v) {
        # Model with only this variable
        formula_only <- as.formula(paste(response_col, "~", v))
        model_only <- maxnet(p = y_true, data = data[, vars], f = formula_only)
        pred_only <- predict(model_only, data[, vars], type = type)
        auc_only <- as.numeric(pROC::auc(y_true, pred_only))
        
        # Model without this variable
        vars_minus_v <- setdiff(vars, v)
        formula_minus <- as.formula(paste(response_col, "~", paste(vars_minus_v, collapse = " + ")))
        model_minus <- maxnet(p = y_true, data = data[, vars], f = formula_minus)
        pred_minus <- predict(model_minus, data[, vars], type = type)
        auc_without <- as.numeric(pROC::auc(y_true, pred_minus))
        
        data.frame(variable = v,
                   auc_only = auc_only,
                   auc_without = auc_without,
                   fold = res$k)
      })
      
      do.call(rbind, fold_res)
      
    } else {
      NULL
    }
  })
  
  # Bind all folds
  jack_all <- dplyr::bind_rows(jack_list)
  
  # Aggregate across folds
  jack_summary <- jack_all %>%
    group_by(variable) %>%
    summarise(mean_auc_only    = mean(auc_only),
              mean_auc_without = mean(auc_without),
              sd_only          = sd(auc_only),
              sd_without       = sd(auc_without),
              n_folds          = n(),
              .groups = "drop") %>%
    arrange(desc(mean_auc_only)) %>%
    mutate(rank = row_number())
  
  return(list(jack_foldwise = jack_all, jack_summary = jack_summary))
}


#Function to compute permutation
permutation_importance <- function(results_list, 
                                   response_col = "presence", 
                                   type = "cloglog", 
                                   n_perm = 10, 
                                   metric = "auc") {
  metric <- match.arg(metric)
  
  # Initialize list
  perm_list <- lapply(results_list, function(res) {
    if(!res$skipped) {
      model <- res$MaxEnt_model
      data  <- res$train_data
      
      vars <- colnames(data)
      vars <- vars[vars != response_col]
      
      # Predict original
      pred_orig <- predict(model, data[, vars], type = type)
      y_true <- data[[response_col]]
      
      # Compute permutation importance
      fold_res <- lapply(vars, function(v) {
        perm_scores <- numeric(n_perm)
        for(i in 1:n_perm) {
          data_perm <- data
          data_perm[[v]] <- sample(data_perm[[v]])
          pred_perm <- predict(model, data_perm[, vars], type = type)
          perm_scores[i] <- switch(metric,
                                   correlation = 1 - cor(pred_orig, pred_perm, use = "complete.obs"),
                                   auc = as.numeric(pROC::auc(y_true, pred_orig)) - as.numeric(pROC::auc(y_true, pred_perm))
          )
        }
        data.frame(variable = v,
                   perm_importance = mean(perm_scores),
                   sd = sd(perm_scores),
                   fold = res$k)
      })
      
      do.call(rbind, fold_res)
    } else {
      NULL
    }
  })
  
  # Bind all folds
  perm_all <- dplyr::bind_rows(perm_list)
  
  # Aggregate across folds
  perm_summary <- perm_all %>%
    group_by(variable) %>%
    summarise(mean_perm = mean(perm_importance),
              sd_perm   = mean(sd),
              n_folds   = n(),
              .groups = "drop") %>%
    arrange(desc(mean_perm)) %>%
    mutate(rank = row_number())
  
  return(list(perm_foldwise = perm_all, perm_summary = perm_summary))
}

#Function to extract responses for plotting
extract_responses <- function(mod, 
                              data_train, 
                              n_points = 500, 
                              id = NULL, 
                              sample_n = NULL) {
  vars <- names(data_train)
  results <- list()
  
  # Optionally sample rows for averaging
  if(!is.null(sample_n) && sample_n < nrow(data_train)) {
    set.seed(42)  # for reproducibility
    data_train <- data_train[sample(nrow(data_train), sample_n), , drop = FALSE]
  }
  
  for(v in vars) {
    # Determine sequence of values for this variable
    if(is.numeric(data_train[[v]])) {
      x_seq <- seq(min(data_train[[v]], na.rm = TRUE),
                   max(data_train[[v]], na.rm = TRUE),
                   length.out = n_points)
    } else {
      x_seq <- sort(unique(data_train[[v]]))
    }
    
    # Numeric variable: predict in batches
    if(is.numeric(data_train[[v]])) {
      tmp_list <- lapply(x_seq, function(xval) {
        tmp <- data_train
        tmp[[v]] <- xval
        tmp
      })
      tmp_all <- do.call(rbind, tmp_list)
      
      pred_all <- predict(mod, newdata = tmp_all, type = "cloglog", clamp = TRUE)
      n_rows <- nrow(data_train)
      y <- colMeans(matrix(pred_all, nrow = n_rows))
      
    } else {
      # Factor variable: loop over levels (usually fast)
      y <- sapply(x_seq, function(xval) {
        tmp <- data_train
        tmp[[v]] <- xval
        mean(predict(mod, newdata = tmp, type = "cloglog", clamp = TRUE))
      })
    }
    
    # Build output
    df <- data.frame(
      variable = v,
      x = x_seq,
      y = y,
      stringsAsFactors = FALSE
    )
    
    if(!is.null(id)) df$id <- id
    results[[v]] <- df
  }
  
  dplyr::bind_rows(results)
}

# 4.3 MaxEnt parallel computing ========================================

## Save predictors final for parallel correct functioning
pred_file_raster <- file.path("Results/Variables_cor", "predictors_final.tif")
writeRaster(predictors_final, pred_file_raster, overwrite=TRUE)

## Loop for sites subsets
for (i in names(list_punts)){ 
  
  #Folder creation for storing results
  out_dir <- file.path("Results/",i,"MaxEnt")
  dir.create(out_dir) #General folder
  dir.create(file.path(out_dir, "Pred_maps")) #Prediction maps
  
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
  
  
  ## Extract the folds in Spatial Block object created in the previous section
  folds <- kfolds$folds_list
  
  ## Parallel loop set up
  ncores <- length(kfolds$folds_list)
  clus <- makeCluster(ncores)
  
  #Load packages inside workers
  clusterEvalQ(clus, {
    library(terra)
    library(ggplot2)
    library(plotROC)
    library(GA)
    library(maxnet)
    library(dplyr)
    library(tidyr)
  })
  
  #Export only simple objects
  to_export <- c("folds", "all_data", "sample_values", "out_dir", "plot_variable_importance")
  clusterExport(clus, varlist = to_export)
  
  #Workers function set up
  worker_kfold_safe <- function(k) {
    
    fold_dir <- file.path(out_dir, as.character(k))
    dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "Pred_maps"), recursive = TRUE, showWarnings = FALSE)
    
    predictors_final <- rast(file.path("Results/Variables_cor", "predictors_final.tif"))
    
    #Train/test split
    trainSet <- unlist(folds[[k]][1])
    testSet  <- unlist(folds[[k]][2])
    
    p <- all_data$Id[trainSet]
    data <- all_data[trainSet, 1:ncol(sample_values)]
    data_test <- all_data[testSet, ]
    
    #GA optimization
    GA_MaxEnt <- ga(
      type = "real-valued",
      fitness = function(params) {
        
        regmult <- params[1]
        fc_index <- max(1, min(round(params[2]), 5))
        fc_comb <- paste(c("l","q","h","p","t")[1:fc_index], collapse="")
        
        # Fit model
        form <- maxnet.formula(p, data, classes = fc_comb)
        mod <- maxnet(p, data, form, regmult=regmult, addsamplestobackground=TRUE)
        
        # serialize
        mod <- unserialize(serialize(mod, NULL))
        
        # AUC on held-out test partition
        preds <- predict(mod, data_test, type="cloglog")
        auc_ME <- as.numeric(pROC::auc(response=data_test$Id, predictor=as.numeric(preds)))
        
        if(is.na(auc_ME)) return(-9999)
        
        # Return AUC directly
        return(auc_ME)
      },
      lower = c(0.5,1),
      upper = c(5.0,5),
      popSize = 20,
      pmutation = 0.4,
      pcrossover = 0.6,
      maxiter = 5,
      elitism = 2,
      seed = 123
    )
    
    best_params <- GA_MaxEnt@solution
    best_regmult <- best_params[1]
    best_fc_num <- max(1, min(round(best_params[2]), 5))
    best_fc <- paste(c("l","q","h","p","t")[1:best_fc_num], collapse="")
    
    # Train final MaxEnt model
    MaxEnt_model <- maxnet(p, data, maxnet.formula(p, data, classes=best_fc),
                           regmult=best_regmult, addsamplestobackground=TRUE)
    
    if(length(MaxEnt_model$betas)==0) return(list(k=k, skipped=TRUE))
    
    # Predict on test set 
    data_test$pred <- predict(MaxEnt_model, data_test, type="cloglog")
    auc_ME <- as.numeric(pROC::auc(response=data_test$Id, predictor=data_test$pred))
    
    # Prediction map creation
    pred_map <- predict(predictors_final, MaxEnt_model, clamp=FALSE, type="cloglog", na.rm=TRUE)
    pred_file <- file.path(out_dir, "Pred_maps", sprintf("pred_map_k%d.tif", k))
    writeRaster(pred_map, pred_file, overwrite=TRUE)
    
    # Response curves 
    resp_file <- file.path(fold_dir, sprintf("Model_response_k%d.tif", k))
    tiff(resp_file, width=6*300, height=8*300, res=300)
    plot(MaxEnt_model, vars=names(MaxEnt_model$samplemeans), common.scale=TRUE, type="cloglog", ylab="Prediction")
    dev.off()
    
    # AUC plot
    auc_file <- file.path(fold_dir, sprintf("AUC_k%d.tif", k))
    ggROC <- ggplot(data_test, aes(m=as.numeric(pred), d=Id)) + geom_roc(n.cuts=0, color='red') + theme_bw()
    tiff(auc_file, width=6*300, height=4*300, res=300)
    print(ggROC + ggtitle(paste("AUC =", round(auc_ME,4))))
    dev.off()
    
    # Return results 
    list(
      k = k,
      auc = auc_ME,
      MaxEnt_model = MaxEnt_model,
      train_data = data,
      test_data = data_test,
      p = p,
      skipped = FALSE
    )
    
  }
  
  ## Collecting Parallel loop results
  results <- parLapplyLB(clus, 1:length(folds), worker_kfold_safe)
  
  # Finishing parallel loop
  stopCluster(clus) 

  save(results, file = file.path(out_dir,paste0("MaxEnt_models.Rdata")))
  
}



# 4.4 Model evaluation =========================================================

## Data frames to store AUC results
auc_df <- data.frame(
  AUC_1 = numeric(length(list_punts)),
  AUC_2 = numeric(length(list_punts)),
  AUC_3 = numeric(length(list_punts)),
  AUC_4 = numeric(length(list_punts)),
  AUC_5 = numeric(length(list_punts)),
  row.names = names(list_punts)
)

##Loop to run over all models
for (i in names(list_punts)){ 
  
##Set folder where to retrieve models
out_dir <- file.path("Results/",i,"MaxEnt")

#Load the models
results <- load(results, file = file.path(out_dir,paste0("MaxEnt_models.Rdata")))

## Creation of mean Prediction map for all folds
pred_maps_list <-  list.files(file.path(out_dir,"Pred_maps"), pattern = "\\.tif[f]?$", full.names = TRUE)
pred_maps_rast <- lapply(pred_maps_list, rast) 
pred_maps_rast <-  do.call(c, pred_maps_rast)
mean_pred_rast <- mean(pred_maps_rast, na.rm = TRUE)

#Save the results
writeRaster(mean_pred_rast, file.path(out_dir,paste0("pred_map_total.tiff")), overwrite=TRUE)

##Response curves for all folds
#Extract values
all_resp <- bind_rows(
  lapply(seq_along(results), function(k) {
    res <- results[[k]]
    mod <- res$MaxEnt_model
    data_train <- res$train_data  
    extract_responses(mod, data_train, id=k)
  })
)

#Plot the results
p <- ggplot() +
  geom_line(data=all_resp, aes(x=x, y=y, group=id, color=factor(id)), alpha=0.7) +
  facet_wrap(~variable, scales="free_x") +
  theme_bw(base_size=14) +
  labs(
    y = "Prediction (cloglog)",
    x = "Variable value",
    color = "Fold"
  ) +
  theme(legend.position="bottom")

ggsave(filename = file.path(out_dir, "Response_Curves.tiff"),
       plot = p,
       width = 12, height = 8, dpi = 300)

##Permutation computation for all folds
perm_results <- permutation_importance(results)

#Plot result
perm_results$perm_summary <- perm_results$perm_summary %>%
  arrange(desc(mean_perm)) %>%
  mutate(variable = factor(variable, levels = variable))

p2 <- ggplot(perm_results$perm_summary, aes(x = variable, y = mean_perm)) +
  geom_col(fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_perm - sd_perm, ymax = mean_perm + sd_perm),
                width = 0.3, color = "black") +
  coord_flip() +  # horizontal bars
  theme_bw(base_size = 14) +
  labs(
    title = "Permutation Importance of Predictors",
    x = "Variable",
    y = "Mean Permutation Importance (± SD)"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(size = 12)
  )

ggsave(file.path(out_dir,"Permutation_Importance.png"), plot = p2,  width = 12, height = 8, dpi = 300)

#save the data
write.csv(perm_results$perm_summary, file = file.path(out_dir,"Permutation_Importance.csv"), )

##Jackknife test for all folds
# occ: vector of 1 (presence) and 0 (background)
# env_data: data.frame or matrix of environmental variables (same rows as occ)
jackknife_results <- jackknife_test(results)

#Plot the results
df <- jack_results$jack_summary %>%
  mutate(variable = factor(variable, levels = variable))

p3 <- ggplot(df, aes(x=variable)) +
  geom_bar(aes(y=only), stat="identity", fill="blue", alpha=0.6) +
  geom_bar(aes(y=without), stat="identity", fill="red", alpha=0.6) +
  ylab("AUC") +
  theme_minimal() +
  ggtitle("Jackknife test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir,"Jackknife.png"), plot = p,  width = 12, height = 8, dpi = 300)

#save the data
write.csv(jackknife_results, file = file.path(out_dir,"Jackknife.csv"), )

##Storing parallel loop results
#AUC
auc_values <- sapply(results[1:length(folds)], function(x) if(x$skipped) NA else x$auc)
auc_df[i, ] <- c(auc_values, rep(NA, ncol(auc_df) - length(auc_values)))

}

## Compute AUC results
auc_df$Mean_AUC <- rowMeans(auc_df[, 1:5], na.rm = TRUE)

#Save results
write.csv(auc_df, file = "Results/AUC_MaxEnt.csv", row.names = TRUE)