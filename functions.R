# adapted from https://github.com/rsetter/opposingclimatevectors

library(terra)
library(data.table)
library(parallel)
library(doParallel)

source("config_dir.R")

#modified to accept climate envelopes rather than ranges
analogcalc <- function(sp_coords, future_rasters, sp_envelope, mask, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  start_time <- Sys.time()
  
  #adjust x to match climate data: -180,180 to 0,360
  sp_coords$x <- ifelse(sp_coords$x < 0, sp_coords$x + 360, sp_coords$x)
  
  #calculate pixel areas
  area_raster <- cellSize(mask,unit='km')
  sp_coords <- as.data.table(sp_coords)
  coords_matrix <- sp_coords[, .(x, y)]
  sp_coords$area <- terra::extract(area_raster, as.matrix(coords_matrix))$area
  sp_coords$sst_future <- terra::extract(future_rasters[[1]], as.matrix(coords_matrix))[[1]]
  sp_coords$ph_future <- terra::extract(future_rasters[[2]], as.matrix(coords_matrix))[[1]]
  
  # Verify future rasters
  if(!all(sapply(future_rasters, inherits, "SpatRaster"))) {
    stop("future_rasters must be a list of two SpatRaster objects (SST and pH)")
  }
  
  # Create unique points database from species coordinates
  cat("ANALOGCALC: Creating points database from species distribution...\n")
  unique_points <- unique(sp_coords[, .(x, y, species_id)])
  setkey(unique_points, x, y)  
  unique_points[, point_id := .I]
  unique_species <- unique(sp_coords$species_id)
  
  # Set number of cores
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  
  # Create wrapped version of mask that can be read in parallel
  cat("Creating memory-mapped objects...\n")
  mask_wrapped <- wrap(mask)
  future_wrapped <- lapply(future_rasters, wrap)
  
  # Calculate optimal chunk size
  total_points <- nrow(unique_points)
  chunk_size <- min(500, max(250, floor(total_points/20)))
  n_chunks <- ceiling(total_points/chunk_size)
  
  cat(sprintf("Processing %d points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("unique_points", "sp_coords","sp_envelope","mask_wrapped","future_wrapped"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  
  # Process chunk of points
  process_chunk <- function(chunk_info) {
    chunk_start <- chunk_info$start
    chunk_end <- chunk_info$end
    
    # Initialize results list for SST and pH
    results_list <- list(
      # SST results
      data.table(
        species_id = integer(), 
        px = numeric(),
        py = numeric(),
        pval = numeric(),
        pvalf = numeric(),
        envelope_min = numeric(),
        envelope_max = numeric(),
        fx = numeric(),
        fy = numeric(),
        fval = numeric(),
        distance = numeric(),
        directionr = numeric(),
        directiondeg = numeric(),
        area = numeric()
      ),
      # pH results
      data.table(
        species_id = integer(), 
        px = numeric(),
        py = numeric(),
        pval = numeric(),
        pvalf = numeric(),
        envelope_min = numeric(),
        envelope_max = numeric(),
        fx = numeric(),
        fy = numeric(),
        fval = numeric(),
        distance = numeric(),
        directionr = numeric(),
        directiondeg = numeric(),
        area = numeric()
      )
    )
    
    # Unwrap mask and future rasters
    mask <- unwrap(mask_wrapped)
    future_rasters <- lapply(future_wrapped, unwrap)
    temp_mask <- mask
    
    # Process each point in the chunk
    for(point_idx in chunk_start:chunk_end) {
      current_point <- unique_points[point_idx]
      
      # Get all species data for this point
      point_data <- sp_coords[x == current_point$x & y == current_point$y]

      # Process each unique species at this point
      unique_species_at_point <- unique(point_data$species_id)
      
      for(sp_id in unique_species_at_point) {
        # Get envelope for current species
        species_envelope <- sp_envelope[sp_envelope$id_no == sp_id,]
        
        if(nrow(species_envelope) == 0) next
        
        # Get data for current species at this point
        species_point_data <- point_data[point_data$species_id == sp_id,][1]
        
      # Process SST and pH
      for(i in 1:2) {  # 1 = SST, 2 = pH
        # Set up current values and envelope bounds
        current_vals <- if(i == 1) {
          list(
            present = species_point_data$sst_mean,
            future = species_point_data$sst_future,
            envelope_min = species_envelope$sst_min,
            envelope_max = species_envelope$sst_max
          )
        } else {
          list(
            present = species_point_data$ph_mean,
            future = species_point_data$ph_future,
            envelope_min = species_envelope$ph_min,
            envelope_max = species_envelope$ph_max
          )
        }
        
          # Extract all future values from raster
          future_vals <- values(future_rasters[[i]])
          future_cells <- which(!is.na(future_vals))
          
          # Find cells within envelope range
          matches <- future_vals[future_cells] <= current_vals$envelope_max & 
            future_vals[future_cells] >= current_vals$envelope_min
          
          if(any(matches)) {
            # Reset temp_mask
            temp_mask[] <- mask[]
            origin_cell <- cellFromXY(temp_mask, 
                                      matrix(c(current_point$x, current_point$y), ncol=2))
            temp_mask[origin_cell] <- 500
            
            dist_raster <- gridDist(temp_mask, target=500, scale=1000)
            
            # Get coordinates of matching cells
            match_cells <- future_cells[matches]
            match_xy <- xyFromCell(future_rasters[[i]], match_cells)
            
            # Calculate distances
            distances <- extract(dist_raster, match_xy)$landseamask
            valid_distances <- !is.na(distances)
            
            if(any(valid_distances)) {
              min_distance <- min(distances[valid_distances])
              closest_idx <- which(distances == min_distance)
              
              for(k in closest_idx) {
                # Calculate direction
                direction <- atan2(match_xy[k,2] - current_point$y,
                                   match_xy[k,1] - current_point$x)
                direction_degrees <- (450 - direction * (180/pi)) %% 360
                
                # Create new row
                new_row <- data.table(
                  species_id = sp_id,
                  px = current_point$x,
                  py = current_point$y,
                  pval = current_vals$present,
                  pvalf = current_vals$future,
                  envelope_min = current_vals$envelope_min,
                  envelope_max = current_vals$envelope_max,
                  fx = match_xy[k,1],
                  fy = match_xy[k,2],
                  fval = future_vals[match_cells[k]],
                  distance = distances[k],
                  directionr = direction,
                  directiondeg = direction_degrees,
                  area = species_point_data$area
                )
                
                results_list[[i]] <- rbindlist(list(results_list[[i]], new_row))
              }
              } else {
                # All distances were NA
                new_row <- data.table(
                  species_id = sp_id,
                  px = current_point$x,
                  py = current_point$y,
                  pval = current_vals$present,
                  pvalf = current_vals$future,
                  envelope_min = current_vals$envelope_min,
                  envelope_max = current_vals$envelope_max,
                  fx = NA_real_,
                  fy = NA_real_,
                  fval = NA_real_,
                  distance = NA_real_,
                  directionr = NA_real_,
                  directiondeg = NA_real_,
                  area = species_point_data$area
                )
                results_list[[i]] <- rbindlist(list(results_list[[i]], new_row))
            }
            } else {
              # No matches within envelope 
              new_row <- data.table(
                species_id = sp_id,
                px = current_point$x,
                py = current_point$y,
                pval = current_vals$present,
                pvalf = current_vals$future,
                envelope_min = current_vals$envelope_min,
                envelope_max = current_vals$envelope_max,
                fx = NA_real_,
                fy = NA_real_,
                fval = NA_real_,
                distance = NA_real_,
                directionr = NA_real_,
                directiondeg = NA_real_,
                area = species_point_data$area
              )
              results_list[[i]] <- rbindlist(list(results_list[[i]], new_row))
        }
      }
      }
      # Garbage collection every 100 points
      if(point_idx %% 100 == 0) {
        rm(dist_raster)
        gc(verbose = FALSE)
      }
    }
    
    return(results_list)
  }
  
  # Process chunks
  chunks <- lapply(1:n_chunks, function(i) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, total_points)
    list(start = start_idx, end = end_idx)
  })
  
  results <- parLapply(cl, chunks, process_chunk)

  stopCluster(cl)
  
  # Combine results
  cat("\nCombining results...\n")
  final_results <- list(
    sst = rbindlist(lapply(results, function(x) x[[1]])),
    ph = rbindlist(lapply(results, function(x) x[[2]]))
  )
  
  end_time <- Sys.time()
  total_time <- as.numeric(end_time - start_time, units="secs")
  cat(sprintf("\nProcessing completed in %0.1f seconds (%0.1f minutes)\n", 
              total_time, total_time/60))

  return(final_results)
}




# Modified analogdifference for envelopes
analogdifference <- function(tempdf, phdf, mask, future_temp, future_ph, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  start_time <- Sys.time()
  
  # Pre-process future values into data.tables for faster lookup
  cat("ANALOGDIFF: Pre-processing future values...\n")
  future_temp_dt <- as.data.table(terra::as.data.frame(future_temp, xy = TRUE))
  setnames(future_temp_dt, c("x", "y", "temp_val"))
  setkey(future_temp_dt, x, y)  # Create index for faster lookups
  
  future_ph_dt <- as.data.table(terra::as.data.frame(future_ph, xy = TRUE))
  setnames(future_ph_dt, c("x", "y", "ph_val"))
  setkey(future_ph_dt, x, y)
  
  # Create database of unique coordinates
  cat("Creating shared points database...\n")
  all_coords <- unique(rbind(
    tempdf[, .(px, py)],
    phdf[, .(px, py)]
  ))
  setnames(all_coords, c("x", "y"))
  all_coords[, point_id := .I]
  
  # Get unique species IDs
  unique_species <- unique(c(tempdf$species_id, phdf$species_id))
  cat(sprintf("Processing %d unique coordinates for %d species...\n", 
              nrow(all_coords), length(unique_species)))
  
  # Set number of cores
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  
  # Create wrapped version of mask for parallel processing
  cat("Creating memory-mapped mask...\n")
  mask_wrapped <- wrap(mask)
  
  # Calculate optimal chunk size
  total_points <- nrow(all_coords)
  chunk_size <- min(500, max(250, floor(total_points/20)))
  n_chunks <- ceiling(total_points/chunk_size)
  
  cat(sprintf("Processing %d unique points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("all_coords", "tempdf", "phdf", "mask_wrapped","future_temp_dt", "future_ph_dt"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  
  # Process chunk of points
  process_chunk <- function(chunk_info) {
    chunk_start <- chunk_info$start
    chunk_end <- chunk_info$end
    
    # Initialize results
    results <- data.table(
      species_id = integer(),  
      px_temp = numeric(),
      py_temp = numeric(),
      pval_temp = numeric(),
      pvalf_temp = numeric(),
      pmin_temp = numeric(),
      pmax_temp = numeric(),
      fx_temp = numeric(),
      fy_temp = numeric(),
      fval_temp = numeric(),
      distance_temp = numeric(),
      directionr_temp = numeric(),
      directiondeg_temp = numeric(),
      px_ph = numeric(),
      py_ph = numeric(),
      pval_ph = numeric(),
      pvalf_ph = numeric(),
      pmin_ph = numeric(),
      pmax_ph = numeric(),
      fx_ph = numeric(),
      fy_ph = numeric(),
      fval_ph = numeric(),
      distance_ph = numeric(),
      directionr_ph = numeric(),
      directiondeg_ph = numeric(),
      area = numeric(),
      analog_diff = numeric(),
      fanalogph_temp = numeric(),
      fanalogphdiff_temp = numeric(),
      fanalogtemp_ph = numeric(),
      fanalogtempdiff_ph = numeric()
    )
    
    # Unwrap mask once at the start
    mask <- unwrap(mask_wrapped)
    temp_mask <- mask
    
    # Process each point in the chunk
    for(point_idx in chunk_start:chunk_end) {
      current_point <- all_coords[point_idx]
      
      # Get all species data for this point
      point_data <- rbind(
        tempdf[px == current_point$x & py == current_point$y],
        phdf[px == current_point$x & py == current_point$y]
      )
      
      # Process each unique species at this point
      unique_species_at_point <- unique(point_data$species_id)
      
      for(sp_id in unique_species_at_point) {
      
      # Get temperature and pH data for current location
      temp_data <- tempdf[px == current_point$x & 
                            py == current_point$y & 
                            species_id == sp_id]
      ph_data <- phdf[px == current_point$x & 
                        py == current_point$y & 
                        species_id == sp_id]
      
      # Skip if either dataset is missing this location
      if(nrow(temp_data) == 0 || nrow(ph_data) == 0) {
        new_row <- create_na_row(current_point, temp_data, ph_data,sp_id)
        results <- rbindlist(list(results, new_row))
        next
      }
      
      # Skip if no future analogs exist
      if(all(is.na(temp_data$fx)) || all(is.na(ph_data$fx))) {
        new_row <- create_na_row(current_point, temp_data, ph_data,sp_id)
        results <- rbindlist(list(results, new_row))
        next
      }
      
      min_distance <- Inf
      best_combo <- NULL
      
      # Calculate distances between temperature and pH analogs
      for(i in seq_len(nrow(temp_data))) {
        if(is.na(temp_data$fx[i])) next
        
        # Reset temp_mask
        temp_mask[] <- mask[]
        tempf_cell <- cellFromXY(temp_mask, 
                                 matrix(c(temp_data$fx[i], temp_data$fy[i]), ncol=2))
        temp_mask[tempf_cell] <- 500
        
        dist_raster <- gridDist(temp_mask, target=500, scale=1000)
        
        for(j in seq_len(nrow(ph_data))) {
          if(is.na(ph_data$fx[j])) next
          
          distance <- extract(dist_raster,matrix(c(ph_data$fx[j], ph_data$fy[j]), ncol=2))$landseamask
          
          if(!is.na(distance) && distance <= min_distance) {
            min_distance <- distance
            
            # Extract future pH at temperature analog location
            fanalogph_temp <- future_ph_dt[.(temp_data$fx[i], temp_data$fy[i]), ph_val]
            
            # Extract future temperature at pH analog location
            fanalogtemp_ph <- future_temp_dt[.(ph_data$fx[j], ph_data$fy[j]), temp_val]
            
            best_combo <- create_result_row(temp_data[i], ph_data[j], distance,fanalogph_temp, fanalogtemp_ph,sp_id)
          }
        }
      }
      
      # Add results
      if(is.null(best_combo)) {
        best_combo <- create_na_row(current_point, temp_data, ph_data,sp_id)
      }
      results <- rbindlist(list(results, best_combo))
      }
      # Garbage collection every 100 points
      if(point_idx %% 100 == 0) {
        rm(dist_raster)
        gc(verbose = FALSE)
      }
    }
    
    return(results)
  }
  
  # Helper function to create NA row
  create_na_row <- function(point, temp_data, ph_data,sp_id) {
    data.table(
      species_id = sp_id,  
      px_temp = point$x,
      py_temp = point$y,
      pval_temp = if(nrow(temp_data) > 0) temp_data$pval[1] else NA_real_,
      pvalf_temp = if(nrow(temp_data) > 0) temp_data$pvalf[1] else NA_real_,
      pmin_temp = if(nrow(temp_data) > 0) temp_data$envelope_min[1] else NA_real_,
      pmax_temp = if(nrow(temp_data) > 0) temp_data$envelope_max[1] else NA_real_,
      fx_temp = if(nrow(temp_data) > 0) temp_data$fx[1] else NA_real_,
      fy_temp = if(nrow(temp_data) > 0) temp_data$fy[1] else NA_real_,
      fval_temp = if(nrow(temp_data) > 0) temp_data$fval[1] else NA_real_,
      distance_temp = if(nrow(temp_data) > 0) temp_data$distance[1] else NA_real_,
      directionr_temp = if(nrow(temp_data) > 0) temp_data$directionr[1] else NA_real_,
      directiondeg_temp = if(nrow(temp_data) > 0) temp_data$directiondeg[1] else NA_real_,
      px_ph = point$x,
      py_ph = point$y,
      pval_ph = if(nrow(ph_data) > 0) ph_data$pval[1] else NA_real_,
      pvalf_ph = if(nrow(ph_data) > 0) ph_data$pvalf[1] else NA_real_,
      pmin_ph = if(nrow(ph_data) > 0) ph_data$envelope_min[1] else NA_real_,
      pmax_ph = if(nrow(ph_data) > 0) ph_data$envelope_max[1] else NA_real_,
      fx_ph = if(nrow(ph_data) > 0) ph_data$fx[1] else NA_real_,
      fy_ph = if(nrow(ph_data) > 0) ph_data$fy[1] else NA_real_,
      fval_ph = if(nrow(ph_data) > 0) ph_data$fval[1] else NA_real_,
      distance_ph = if(nrow(ph_data) > 0) ph_data$distance[1] else NA_real_,
      directionr_ph = if(nrow(ph_data) > 0) ph_data$directionr[1] else NA_real_,
      directiondeg_ph = if(nrow(ph_data) > 0) ph_data$directiondeg[1] else NA_real_,
      area = if(nrow(temp_data) > 0) temp_data$area[1] else NA_real_,
      analog_diff = NA_real_,
      fanalogph_temp = NA_real_,
      fanalogphdiff_temp = NA_real_,
      fanalogtemp_ph = NA_real_,
      fanalogtempdiff_ph = NA_real_
    )
  }
  
  # Helper function to create result row
  create_result_row <- function(temp_data, ph_data, distance, fanalogph_temp, fanalogtemp_ph,sp_id) {
    data.table(
      species_id = sp_id, 
      px_temp = temp_data$px,
      py_temp = temp_data$py,
      pval_temp = temp_data$pval,
      pvalf_temp = temp_data$pvalf,
      pmin_temp = temp_data$envelope_min,
      pmax_temp = temp_data$envelope_max,
      fx_temp = temp_data$fx,
      fy_temp = temp_data$fy,
      fval_temp = temp_data$fval,
      distance_temp = temp_data$distance,
      directionr_temp = temp_data$directionr,
      directiondeg_temp = temp_data$directiondeg,
      px_ph = ph_data$px,
      py_ph = ph_data$py,
      pval_ph = ph_data$pval,
      pvalf_ph = ph_data$pvalf,
      pmin_ph = ph_data$envelope_min,
      pmax_ph = ph_data$envelope_max,
      fx_ph = ph_data$fx,
      fy_ph = ph_data$fy,
      fval_ph = ph_data$fval,
      distance_ph = ph_data$distance,
      directionr_ph = ph_data$directionr,
      directiondeg_ph = ph_data$directiondeg,
      area = temp_data$area,
      analog_diff = distance,
      fanalogph_temp = fanalogph_temp,
      fanalogphdiff_temp = if(length(fanalogph_temp) > 0) fanalogph_temp - ph_data$pval else NA_real_,
      fanalogtemp_ph = fanalogtemp_ph,
      fanalogtempdiff_ph = if(length(fanalogtemp_ph) > 0) fanalogtemp_ph - temp_data$pval else NA_real_
    )
  }
  
  # Prepare chunks
  chunks <- lapply(1:n_chunks, function(i) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, total_points)
    list(start = start_idx, end = end_idx)
  })
  
  # Process chunks in parallel
  results <- parLapply(cl, chunks, process_chunk)
  end_time <- Sys.time()
  
  stopCluster(cl)
  
  # Combine results
  cat("\nCombining results...\n")
  final_results <- rbindlist(results)
  
  # Verify results
  expected_rows <- length(unique_species) * nrow(all_coords)
  actual_rows <- nrow(final_results)
  cat(sprintf("Expected %d rows, got %d rows\n", 
              expected_rows, actual_rows))
  
  total_time <- as.numeric(end_time - start_time, units="secs")
  cat(sprintf("\nProcessing completed in %0.1f seconds (%0.1f minutes)\n", 
              total_time, total_time/60))
  
  return(final_results)
}



# Modified dualanalog for envelopes
dualanalog <- function(sp_coords, sp_envelope, future_temp, future_ph, mask, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  cat("DUALCALC: Initializing data structures...\n")
  start_time <- Sys.time()
  
  terra::terraOptions(memfrac=0.5)
  
  # Validate inputs
  if (!all(sapply(list(future_temp, future_ph), 
                  function(x) inherits(x, "SpatRaster")))) {
    stop("future_temp and future_ph must be SpatRaster objects")
  }
  
  # Convert sp_coords to data.table if it isn't already
  if (!is.data.table(sp_coords)) {
    sp_coords <- as.data.table(sp_coords)
  }
  if (!is.data.table(sp_envelope)) {
    sp_envelope <- as.data.table(sp_envelope)
  }
  
  # Adjust x coordinates if needed (-180,180 to 0,360)
  sp_coords[, x := ifelse(x < 0, x + 360, x)]
  
  # Extract future conditions for species locations
  sp_coords[, ':='(
    sst_future = terra::extract(future_temp, cbind(x, y))[[1]],
    ph_future = terra::extract(future_ph, cbind(x, y))[[1]]
  )]
  
  # Create all_data structure
  sp_coords_subset <- sp_coords[, .(
    species_id,
    x, y,
    sst_mean,
    ph_mean,
    sst_future,
    ph_future
  )]
  
  # Merge with envelope data
  all_data <- merge(
    sp_coords_subset,
    sp_envelope[, .(
      species_id = id_no,
      temp_min = sst_min,
      temp_max = sst_max,
      ph_min = ph_min,
      ph_max = ph_max
    )],
    by = "species_id")
  
  # Add IDs
  all_data[, id := .I]
  
  # Set keys for faster operations
  setkey(all_data, x, y)
  
  # Create all_coords from all_data
  all_coords <- all_data[, .(x, y)]
  
  # Set number of cores
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  
  # Create wrapped version of mask for parallel processing
  cat("Creating memory-mapped mask...\n")
  mask_wrapped <- wrap(mask)
  tos_wrapped <- wrap(future_temp)
  ph_wrapped <- wrap(future_ph)
  
  # Calculate optimal chunk size
  total_points <- nrow(all_coords)
  chunk_size <- min(500, max(250, floor(total_points/20)))
  n_chunks <- ceiling(total_points/chunk_size)
  
  cat(sprintf("Processing %d unique points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("all_coords", "all_data", "mask_wrapped","tos_wrapped","ph_wrapped"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  
  # Helper function to create NA row
  create_na_row <- function(data) {
    data.table(
      species_id = data$species_id,
      px = data$x,
      py = data$y,
      pval_temp = data$sst_mean,
      pvalf_temp = data$sst_future,
      pmin_temp = data$temp_min,
      pmax_temp = data$temp_max,
      pval_ph = data$ph_mean,
      pvalf_ph = data$ph_future,
      pmin_ph = data$ph_min,
      pmax_ph = data$ph_max,
      fx = NA_real_,
      fy = NA_real_,
      fval_temp = NA_real_,
      fval_ph = NA_real_,
      distance = NA_real_,
      directionr = NA_real_,
      directiondeg = NA_real_
    )
  }
  
  # Process chunk of points
  process_chunk <- function(chunk_info) {
    chunk_start <- chunk_info$start
    chunk_end <- chunk_info$end
    
    # Initialize results
    results <- data.table(
      species_id = integer(),
      px = numeric(),
      py = numeric(),
      pval_temp = numeric(),
      pvalf_temp = numeric(),
      pmin_temp = numeric(),
      pmax_temp = numeric(),
      pval_ph = numeric(),
      pvalf_ph = numeric(),
      pmin_ph = numeric(),
      pmax_ph = numeric(),
      fx = numeric(),
      fy = numeric(),
      fval_temp = numeric(),
      fval_ph = numeric(),
      distance = numeric(),
      directionr = numeric(),
      directiondeg = numeric()
    )
    
    # Unwrap mask once at the start
    mask <- unwrap(mask_wrapped)
    temp_mask <- mask
    future_temp <- unwrap(tos_wrapped)
    future_ph <- unwrap(ph_wrapped)
                
    # Process each point in the chunk
    for(point_idx in chunk_start:chunk_end) {
      # Get data for current location
      current_data <- all_data[point_idx]
      
      if(nrow(current_data) == 0) {
        results <- rbind(results, create_na_row(current_data))
        next
      }
      
      # Create raster masks for valid regions
      temp_valid <- future_temp <= current_data$temp_max & future_temp >= current_data$temp_min
      ph_valid <- future_ph <= current_data$ph_max & future_ph >= current_data$ph_min
      
      # Combine conditions
      valid_cells <- which(!is.na(temp_valid[]) & temp_valid[] & 
                             !is.na(ph_valid[]) & ph_valid[] & 
                             !is.na(mask[]))
      
      if(length(valid_cells) == 0) {
        results <- rbind(results, create_na_row(current_data))
        next
      }
      
      # Get coordinates of valid cells
      valid_xy <- xyFromCell(future_temp, valid_cells)
      
      # Calculate distances for valid points
      temp_mask <- mask
      origin_cell <- cellFromXY(temp_mask, 
                                matrix(c(current_data$x, current_data$y), ncol=2))
      temp_mask[origin_cell] <- 500
      
      dist_raster <- gridDist(temp_mask, target=500, scale=1000)
      
      distances <- extract(dist_raster, valid_xy)$landseamask
      valid_coords <- data.table(
        x = valid_xy[,1],
        y = valid_xy[,2],
        distance = distances
      )
      valid_coords <- valid_coords[!is.na(distance)]
      
      if(nrow(valid_coords) > 0) {
        closest_match <- valid_coords[which.min(distance)]
        
        # Calculate direction
        direction <- atan2(
          closest_match$y - current_data$y,
          closest_match$x - current_data$x
        )
        direction_degrees <- (450 - direction * (180/pi)) %% 360
        
        closest_temp <- terra::extract(future_temp, 
                                       matrix(c(closest_match$x, closest_match$y), ncol=2))[[1]]
        closest_ph <- terra::extract(future_ph, 
                                     matrix(c(closest_match$x, closest_match$y), ncol=2))[[1]]
        
        # Add result
        results <- rbind(results, data.table(
          species_id = current_data$species_id,
          px = current_data$x,
          py = current_data$y,
          pval_temp = current_data$sst_mean,
          pvalf_temp = current_data$sst_future,
          pmin_temp = current_data$temp_min,
          pmax_temp = current_data$temp_max,
          pval_ph = current_data$ph_mean,
          pvalf_ph = current_data$ph_future,
          pmin_ph = current_data$ph_min,
          pmax_ph = current_data$ph_max,
          fx = closest_match$x,
          fy = closest_match$y,
          fval_temp = closest_temp,
          fval_ph = closest_ph,
          distance = closest_match$distance,
          directionr = direction,
          directiondeg = direction_degrees
        ))
      } else {
        results <- rbind(results, create_na_row(current_data))
      }
      
      # Garbage collection every 100 points
      if(point_idx %% 100 == 0) {
        rm(dist_raster)
        gc(verbose = FALSE)
      }
    }
    
    return(results)
  }
  
  # Prepare chunks
  chunks <- lapply(1:n_chunks, function(i) {
    start_idx <- (i-1) * chunk_size + 1
    end_idx <- min(i * chunk_size, total_points)
    list(start = start_idx, end = end_idx)
  })
  
  # Process chunks in parallel
  results <- parLapply(cl, chunks, process_chunk)
  end_time <- Sys.time()
  
  stopCluster(cl)
  
  # Combine results
  cat("\nCombining results...\n")
  final_results <- rbindlist(results)
  
  # Verify results
  expected_rows <- nrow(all_coords)
  actual_rows <- nrow(final_results)
  cat(sprintf("Expected %d rows, got %d rows\n", 
              expected_rows, actual_rows))
  
  total_time <- as.numeric(end_time - start_time, units="secs")
  cat(sprintf("\nProcessing completed in %0.1f seconds (%0.1f minutes)\n", 
              total_time, total_time/60))
  
  return(final_results)
}



# modified to run per species with envelopes
run_analog_analysis <- function(species_list, species_envelopes, 
                                future_temp_list, future_ph_list, 
                                time_periods, base_period="hist", mask, output_dir,
                                scenario, model, ncores=NULL,
                                batch_size=50) {
  require(terra)
  require(data.table)
  
  # Validate inputs
  if(length(future_temp_list) != length(future_ph_list) || 
     length(future_temp_list) != length(time_periods)) {
    stop("Number of future rasters must match number of time periods")
  }
  
  # Validate base_period
  if(!base_period %in% c("hist", "1985")) {
    stop("base_period must be either 'hist' or '1985'")
  }
  
  # Create temporary directory to store intermediate progress
  temp_dir <- file.path(output_dir, "temp_batches")
  if(!dir.exists(temp_dir)) {
    dir.create(temp_dir)
  }
  
  # Convert inputs to data.tables if not already
  species_dt <- as.data.table(species_list)
  envelopes_dt <- as.data.table(species_envelopes)
  
  # Get unique species IDs
  species_ids <- unique(species_dt$species_id)
  n_species <- length(species_ids)
  
  # Calculate number of batches
  n_batches <- ceiling(n_species / batch_size)
  
  cat(sprintf("Processing %d species in %d batches of %d species each\n", 
              n_species, n_batches, batch_size))
  
  # Process each time period
  for(t in seq_along(time_periods)) {
    end_period <- time_periods[t]
    future_temp <- future_temp_list[[t]]
    future_ph <- future_ph_list[[t]]
    
    cat(sprintf("\nProcessing time period: %s to %s\n", base_period, end_period))
    
    # Initialize lists to store results for this time period
    period_results <- list(
      dual = vector("list", n_batches),
      tos = vector("list", n_batches),
      ph = vector("list", n_batches),
      diff = vector("list", n_batches)
    )
    
    # Process species in batches
    for(b in 1:n_batches) {
      batch_start <- (b-1) * batch_size + 1
      batch_end <- min(b * batch_size, n_species)
      batch_species <- species_ids[batch_start:batch_end]
      
      # Create batch-specific temp filenames
      batch_files <- list(
        dual = file.path(temp_dir, sprintf("batch_%d_dualanalog_%s_%s_%s-%s.csv",
                                           b, model, scenario, base_period, end_period)),
        tos = file.path(temp_dir, sprintf("batch_%d_analog_tos_%s_%s_%s-%s.csv",
                                          b, model, scenario, base_period, end_period)),
        ph = file.path(temp_dir, sprintf("batch_%d_analog_ph_%s_%s_%s-%s.csv",
                                         b, model, scenario, base_period, end_period)),
        diff = file.path(temp_dir, sprintf("batch_%d_analogdiff_%s_%s_%s-%s.csv",
                                           b, model, scenario, base_period, end_period))
      )
      
      cat(sprintf("\nProcessing batch %d/%d (species %d to %d)\n", 
                  b, n_batches, batch_start, batch_end))
      
      # Get data for current batch of species
      batch_coords <- species_dt[species_id %in% batch_species]
      batch_envelopes <- envelopes_dt[id_no %in% batch_species]
      
      # Run analyses for this batch
      # run dual analog
      dual_results <- dualanalog(
        sp_coords = batch_coords,
        sp_envelope = batch_envelopes,
        future_temp = future_temp,
        future_ph = future_ph,
        mask = mask,
        ncores = ncores
      )
      
      # Save batch dual results 
      fwrite(dual_results, batch_files$dual, row.names = FALSE)
      period_results$dual[[b]] <- dual_results
      
      # run analog calc which returns a list of sst and ph results
      tosph_results <- analogcalc(
        sp_coords = batch_coords,
        future_rasters = list(future_temp, future_ph),
        sp_envelope = batch_envelopes,
        mask = mask,
        ncores = ncores
      )
      
      # Save batch tosph results 
      fwrite(tosph_results$sst, batch_files$tos, row.names = FALSE)
      fwrite(tosph_results$ph, batch_files$ph, row.names = FALSE)
      period_results$tos[[b]] <- tosph_results$sst
      period_results$ph[[b]] <- tosph_results$ph
      
      # run analog difference
      diff_results <- analogdifference(
        tempdf = tosph_results$sst,
        phdf = tosph_results$ph,
        mask = mask,
        future_temp = future_temp,
        future_ph = future_ph,
        ncores = ncores
      )
      
      # Save batch diff results 
      fwrite(diff_results, batch_files$diff, row.names = FALSE)
      period_results$diff[[b]] <- diff_results
      
      # Store results for this batch
      period_results$dual[[b]] <- dual_results
      period_results$tos[[b]] <- tosph_results$sst
      period_results$ph[[b]] <- tosph_results$ph
      period_results$diff[[b]] <- diff_results
      
      # Garbage collection
      gc()
    }
    
    # Combine all batches for this time period
    cat("\nCombining batch results...\n")
    combined_results <- list(
      dual = rbindlist(period_results$dual),
      tos = rbindlist(period_results$tos),
      ph = rbindlist(period_results$ph),
      diff = rbindlist(period_results$diff)
    )
    
    # Save combined results for this time period
    cat("Saving results...\n")
    output_files <- list(
      dual = file.path(output_dir, 
                       sprintf("dualanalog_%s_%s_%s-%s.csv",
                               model, scenario, base_period, end_period)),
      tos = file.path(output_dir,
                      sprintf("analog_tos_%s_%s_%s-%s.csv",
                              model, scenario, base_period, end_period)),
      ph = file.path(output_dir,
                     sprintf("analog_ph_%s_%s_%s-%s.csv",
                             model, scenario, base_period, end_period)),
      diff = file.path(output_dir,
                       sprintf("analogdiff_%s_%s_%s-%s.csv",
                               model, scenario, base_period, end_period))
    )
    
    # Save all results
    for(type in names(combined_results)) {
      fwrite(combined_results[[type]], 
             output_files[[type]], 
             row.names = FALSE)
      
      # Print number of rows saved
      cat(sprintf("Saved %d rows to %s\n", 
                  nrow(combined_results[[type]]),
                  basename(output_files[[type]])))
    }
    
    # Clean up this time period's results
    rm(combined_results, period_results)
    gc()
    
    # Clean up temp batch files after successful completion of time period
    unlink(list.files(temp_dir, full.names = TRUE))
  }
}




#notification that code has completed
notifyitsdone <- function(api_token,project = "rcode",channel = "rproject",event = "Analysis Complete",
                          description = "Your R code has finished running",icon = "✅"){
  library(RCurl)
  
  # set up notification system
  headers = c(
    "Content-Type" = "application/json",
    "Authorization" = api_token #insert personal API token (generated when sign up)
  )

  params = "{
  \"project\": \"rcode\",
  \"channel\": \"rproject\",
  \"event\": \"Analysis Complete\",
  \"description\": \"Your R code has finished running\",
  \"icon\": \"🔥\",
  \"notify\": true
}"
  
  # Replace each value while preserving format
  params = gsub('"rproject"', paste0('"', channel, '"'), params, fixed=TRUE)
  params = gsub('"Analysis Complete"', paste0('"', event, '"'), params, fixed=TRUE)
  params = gsub('"Your R code has finished running"', paste0('"', description, '"'), params, fixed=TRUE)
  params = gsub('"🔥"', paste0('"', icon, '"'), params, fixed=TRUE)
  
  #code below sends the notificaiton
  res <- postForm("https://api.logsnag.com/v1/log", .opts=list(postfields = params, httpheader = headers, followlocation = TRUE), style = "httppost") 
  cat("notification sent to phone\n")
}





