##########################################################################
# This script is Point Process models for The settlement pattern         #
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

# 1 Load folders of the subsets ################################################

# Define your directory path
path <- "Results"

# List all items in the directories
folders <- list.dirs(path, recursive = FALSE, full.names = TRUE)

# Filter only those whose names start with "v"
folders <- folders[!grepl("^v", basename(folders), ignore.case = TRUE)]

#Function to conver raster covariates to im
as.im.SpatRaster1 <- function(X) {
  X <- X[[1]]  # Take first layer
  
  if (terra::ncell(X) == 0 || all(is.na(X[]))) {
    stop("Raster has no data or only NA values.")
  }
  
  rs <- terra::res(X)
  e <- as.vector(terra::ext(X))
  
  mat <- as.matrix(X, wide = TRUE)
  if (is.null(mat) || ncol(mat) == 0) {
    stop("Raster matrix is empty or invalid.")
  }
  
  mat <- mat[nrow(mat):1, ]  # Flip matrix for spatstat orientation
  
  is_factor <- terra::is.factor(X)
  
  if (is_factor) {
    levs <- terra::levels(X)[[1]]
    
    if (is.null(levs) || !("ID" %in% colnames(levs))) {
      stop("Factor raster has invalid or missing level definitions.")
    }
    
    labels_column <- if ("category" %in% colnames(levs)) {
      "category"
    } else {
      colnames(levs)[2]
    }
    
    labels <- as.character(levs[[labels_column]])
    mat <- factor(mat, levels = levs$ID, labels = labels)
    
    # ✅ Restore 2D matrix structure for factor
    dim(mat) <- dim(X)[1:2]
    
    im_type <- "factor"
  } else {
    im_type <- "real"
  }
  
  out <- list(
    v = mat,
    dim = dim(X)[1:2],
    xrange = e[1:2],
    yrange = e[3:4],
    xstep = rs[1],
    ystep = rs[2],
    xcol = e[1] + (1:ncol(mat)) * rs[1] - 0.5 * rs[1],
    yrow = e[4] - (nrow(mat):1) * rs[2] + 0.5 * rs[2],
    type = im_type,
    units = list(singular = terra::units(X), plural = terra::units(X), multiplier = 1)
  )
  
  attr(out$units, "class") <- "unitname"
  class(out) <- "im"
  return(out)
} #SpatRaster to im function

# 2 PPMs with euclidean distances and pred map #################################

##Loop for computing PPMs for every folder
for(f in 1:length(folders)) {
  
  #Set directory
  folder <- folders[[f]]
  folder_path <- file.path(folder,"MaxEnt")
  
  #Load the rasters
  file_pred   <- file.path(folder_path, "pred_map_total.tiff")
  
  if (!file.exists(file_pred)) {
    message("Skipping ", folder, ": pred_map_total.tiff not found.")
    next
  }
  
  #Convert raster to im
  raster_pred <- rast(file.path(folder_path, "pred_map_total.tiff"))
  pred_im <- lapply(raster_pred, as.im.SpatRaster1) 
  names(pred_im) <- "Prediction_map"
  
  #Load the points to ppp
  locs <- st_read(file.path(folder,paste0(basename(folder),".shp")))
  data_jac <- spatstat.geom::as.ppp(locs[3:4])
  spatstat.geom::marks(data_jac) <- NULL
  
  ## Euclidean K function 
  K <- Kest(data_jac)
  EK <- envelope(data_jac, Kest, verbose = FALSE, correction = "best")
  
  #Plot the results
  pdf(file.path(folder_path, "K_base.pdf"))
  plot(EK, col = "red",
       main = paste("K function with euclidean distances of", basename(folder)),
       cex.main = 0.9)
  dev.off()
  
  ## Fit the PPM with pred map
  ppm_model <- ppm(data_jac ~ Prediction_map,
                   covariates = pred_im)
  
  #Evaluate Residuals from envelope in parallel
  nsim <- 99  # number of simulations
  ncores <- 4
  sim_per_core <- ceiling(nsim / ncores) # split simulations across cores
  
  #Create cluster
  clus <- makeCluster(ncores)
  
  #Load required packages on each worker
  clusterEvalQ(clus, {
    library(spatstat.explore)
    library(spatstat.model)
  })
  
  #Export ppm_model to workers
  clusterExport(clus, varlist = c("ppm_model","sim_per_core"))
  
  #Function for workers in parallel
  worker_envelope <- function(dummy) {
    envelope(ppm_model, fun = Kres, nsim = sim_per_core, nrank = 1,
             correction = "best", simulate = NULL)
  }
  
  #Run envelopes in parallel
  envelope_chunks <- parLapply(clus, 1:ncores, worker_envelope)
  
  stopCluster(clus) #Stop cluster
  
  ##Collect results of parallel processing
  #Combine results into a single envelope object
  simpatterns_list <- lapply(envelope_chunks, function(x) attr(x, "simpatterns"))
  observed_list <- lapply(envelope_chunks, function(x) x$obs)
  theoretical_list <- lapply(envelope_chunks, function(x) x$theo)
  lower_envelope_list <- lapply(envelope_chunks, function(x) x$lo)
  upper_envelope_list <- lapply(envelope_chunks, function(x) x$hi)
  
  #Combine the individual envelope components
  combined_simpatterns <- do.call(c, simpatterns_list)
  combined_observed <- unlist(observed_list)
  combined_theoretical <- unlist(theoretical_list)
  combined_lower_envelope <- unlist(lower_envelope_list)
  combined_upper_envelope <- unlist(upper_envelope_list)
  
  #Create a new combined envelope object using the combined data
  if (all(combined_simpatterns == 0) && all(combined_theoretical == 0)) {
    combined_envelope <- data.frame(
      r = envelope_chunks[[1]]$r,  # Use 'r' from the first chunk
      obs = combined_observed,
      lo = combined_lower_envelope,
      hi = combined_upper_envelope
    )
  } else {
    # If not empty, combine all components
    combined_envelope <- data.frame(
      r = envelope_chunks[[1]]$r,  # Use 'r' from the first chunk
      obs = combined_observed,
      theo = combined_theoretical,
      lo = combined_lower_envelope,
      hi = combined_upper_envelope
    )
  }
  
  # Add necessary attributes to the new envelope object
  attr(combined_envelope, "class") <- c("envelope", "fv", "data.frame")
  attr(combined_envelope, "argu") <- "r"
  attr(combined_envelope, "valu") <- "obs"
  attr(combined_envelope, "ylab") <- expression(bold(R) ~ hat(K)(r))
  attr(combined_envelope, "fmla") <- ". ~ r"
  attr(combined_envelope, "labl") <- c("r", "hat(%s)[obs](r)", "%s[theo](r)", "hat(%s)[lo](r)", "hat(%s)[hi](r)")
  attr(combined_envelope, "einfo") <- envelope_chunks[[1]]$einfo  # You can adjust this if necessary
  attr(combined_envelope, "desc") <- "K-function Envelope"
  theo_values <- combined_envelope$theo
  
  #Plot the results
  pdf(file.path(folder_path, "Kres_pred_map.pdf"))
  plot(combined_envelope, col = "red", lwd = 1,
       main = paste(basename(folder), "goodness of fit (Kres +", as.character(nsim),"simulations of envelope)"),
       cex.main = 0.9, ylab = expression(bold(R) ~ hat(K)(r)), xlab = "r (distance)", legend = FALSE)
  legend("topleft",
         legend = c(expression(hat(K)[obs](r)),
                    expression(hat(K)[lo](r)),
                    expression(hat(K)[hi](r))),
         col = c("red", "gray", "gray"),
         lty = c(2, 1, 1), lwd = 2, cex = 1,
         bty = "o", box.lwd = 1, box.col = "black")
  dev.off()
  
}

# 3 Clustered PPMs with pred map ###############################################

## Lopp for all three cluster models
for (model in c("Th", "Ma", "Ca")) {
  
  #Assign desired model
  KPPM_model <- model
  
  #Loop for all folders
  for(f in 1:length(folders)) {
    
    #Set directory
    folder <- folders[[f]]
    folder_path <- file.path(folder,"MaxEnt")
    
    #Load the rasters
    file_pred   <- file.path(folder_path, "pred_map_total.tiff")
    
    if (!file.exists(file_pred)) {
      message("Skipping ", folder, ": pred_map_total.tiff not found.")
      next
    }
    
    #Convert raster to im
    raster_pred <- rast(file.path(folder_path, "pred_map_total.tiff"))
    pred_im <- lapply(raster_pred, as.im.SpatRaster1) 
    names(pred_im) <- "Prediction_map"
    
    #Load the points to ppp
    locs <- st_read(file.path(folder,paste0(basename(folder),".shp")))
    data_jac <- spatstat.geom::as.ppp(locs[3:4])
    spatstat.geom::marks(data_jac) <- NULL
    
    ##Euclidean KPPM (Pred_map) with Thomas, Matérn and Cauchy
    if(KPPM_model == "Th") {
      ppm_model <- kppm(X = data_jac, trend = ~ Prediction_map, covariates = pred_im, clusters = "Thomas")
    } else if (KPPM_model == "Ma") {
      ppm_model <- kppm(X = data_jac, trend = ~ Prediction_map, covariates = pred_im, clusters = "MatClus")
    } else if (KPPM_model == "Ca") {
      ppm_model <- kppm(X = data_jac, trend = ~ Prediction_map, covariates = pred_im, clusters = "Cauchy")
    }
    
    #Evaluate KPPM with Kres envelope
    nsim <- 99  # number of simulations
    ncores <- 4
    sim_per_core <- ceiling(nsim / ncores) # split simulations across cores
    
    #Create cluster
    clus <- makeCluster(ncores)
    
    #Load required packages on each worker
    clusterEvalQ(clus, {
      library(spatstat.explore)
      library(spatstat.model)
    })
    
    #Export ppm_model to workers
    clusterExport(clus, varlist = c("ppm_model","sim_per_core"))
    
    #Function for workers in parallel
    worker_envelope <- function(dummy) {
      envelope.kppm(ppm_model, fun = Kres, nsim = sim_per_core)
    }
    
    # Run envelopes in parallel
    envelope_chunks <- parLapply(clus, 1:ncores, worker_envelope)
    
    stopCluster(clus) #Stop cluster
    
    ##Collect results of parallel processing 
    # Combine results into a single envelope object
    simpatterns_list <- lapply(envelope_chunks, function(x) attr(x, "simpatterns"))
    observed_list <- lapply(envelope_chunks, function(x) x$obs)
    theoretical_list <- lapply(envelope_chunks, function(x) x$theo)
    lower_envelope_list <- lapply(envelope_chunks, function(x) x$lo)
    upper_envelope_list <- lapply(envelope_chunks, function(x) x$hi)
    
    # Combine the individual envelope components
    combined_simpatterns <- do.call(c, simpatterns_list)
    combined_observed <- unlist(observed_list)
    combined_theoretical <- unlist(theoretical_list)
    combined_lower_envelope <- unlist(lower_envelope_list)
    combined_upper_envelope <- unlist(upper_envelope_list)
    
    # Create a new combined envelope object using the combined data
    if (all(combined_simpatterns == 0) && all(combined_theoretical == 0)) {
      combined_envelope <- data.frame(
        r = envelope_chunks[[1]]$r,  # Use 'r' from the first chunk
        obs = combined_observed,
        lo = combined_lower_envelope,
        hi = combined_upper_envelope
      )
    } else {
      # If not empty, combine all components
      combined_envelope <- data.frame(
        r = envelope_chunks[[1]]$r,  # Use 'r' from the first chunk
        obs = combined_observed,
        theo = combined_theoretical,
        lo = combined_lower_envelope,
        hi = combined_upper_envelope
      )
    }
    
    # Add necessary attributes to the new envelope object
    attr(combined_envelope, "class") <- c("envelope", "fv", "data.frame")
    attr(combined_envelope, "argu") <- "r"
    attr(combined_envelope, "valu") <- "obs"
    attr(combined_envelope, "ylab") <- expression(bold(R) ~ hat(K)(r))
    attr(combined_envelope, "fmla") <- ". ~ r"
    attr(combined_envelope, "labl") <- c("r", "hat(%s)[obs](r)", "%s[theo](r)", "hat(%s)[lo](r)", "hat(%s)[hi](r)")
    attr(combined_envelope, "einfo") <- envelope_chunks[[1]]$einfo  # You can adjust this if necessary
    attr(combined_envelope, "desc") <- "K-function Envelope"
    theo_values <- combined_envelope$theo
    
    #plot the results
    pdf(file.path(folder_path, paste0("KPPM_",as.character(KPPM_model),".pdf")))
    par(mfrow = c(1, 1), mar = c(3, 3, 3, 3))
    plot(combined_envelope, 
         col = "red",  # Color for the observed K-function
         lwd = 1,      # Line width
         main = paste(basename(folder), "goodness of fit in",as.character(KPPM_model),"(Kres +", as.character(nsim),"simulations of envelope)"),  # Main title
         cex.main = 0.9,  # Font size of the title
         ylab = expression(bold(R) ~ hat(K)(r)),  # Y-axis label
         xlab = "r (distance)",  # X-axis label
         legend = FALSE  # Disable automatic legend
    )
    legend("topleft", 
           legend = c(expression(hat(K)[obs](r)), 
                      expression(hat(K)[lo](r)), 
                      expression(hat(K)[hi](r))),  # Custom labels for the legend
           col = c("red", "gray", "gray"),  # Line colors
           lty = c(2, 1, 1),  # Line types (solid and dashed)
           lwd = 2,  # Line width
           cex = 1,  # Adjust text size
           bty = "o",
           box.lwd = 1,  # Border line width
           box.col = "black"  # Border color
    )
    
    dev.off()
    
  }
  
}
