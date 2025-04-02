# adapted from https://github.com/rsetter/opposingclimatevectors

library(terra)
library(data.table)
library(parallel)
library(doParallel)

source("config_dir.R")
scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
years <- c(2020,2030,2040,2050,2060,2070,2080,2090,2100)


#analog functions. adapted from climate vectors project

#modified to accept climate envelopes rather than ranges
analogcalc <- function(sp_coords, future_temp, future_ph, sp_envelope, mask, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  start_time <- Sys.time()
  
  #increase memory allocation
  terra::terraOptions(memfrac=0.8,progress=10)

  #adjust x to match climate data: -180,180 to 0,360
  sp_coords$x <- ifelse(sp_coords$x < 0, sp_coords$x + 360, sp_coords$x)
  
  #calculate pixel areas
  sp_coords <- as.data.table(sp_coords)
  sp_envelope <- as.data.table(sp_envelope)
  coords_matrix <- sp_coords[, .(x, y)]
  area_raster <- cellSize(mask,unit='km')
  sp_coords$area <- terra::extract(area_raster, as.matrix(coords_matrix))$area
  setkey(sp_coords, x, y)
  setkey(sp_envelope, id_no)

  # Create unique points database from species coordinates
  cat("ANALOGCALC: Creating points database from species distribution...\n")
  unique_points <- unique(sp_coords[, .(x, y)])
  setkey(unique_points, x, y)  
  unique_points[, point_id := .I]
  unique_species <- unique(sp_coords$species_id)
  
  #set up results columns
  result_columns <- c(
    "species_id", "px", "py", "pval", "pvalf", "envelope_min", "envelope_max", 
    "fx", "fy", "fval", "distance", "directionr", "directiondeg", "area"
  )
  
  # Set number of cores and chunks
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  total_points <- nrow(unique_points)
  points_per_chunk <- min(5000,max(250, floor(total_points/(4 * ncores))))
  n_chunks <- ceiling(total_points/points_per_chunk)

  # Create wrapped version of mask that can be read in parallel
  mask_wrapped <- wrap(mask)
  future_temp_wrapped <- wrap(future_temp)
  future_ph_wrapped <- wrap(future_ph)
  
  cat(sprintf("Processing %d points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("unique_points", "sp_coords","sp_envelope","mask_wrapped","future_temp_wrapped", "future_ph_wrapped"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  
  # Helper function to create result row
  create_result_row <- function(point, sp_id, current_vals, area, future_xy = NULL, future_val = NULL, distance = NULL, direction = NULL) {
    data.table(
      species_id = sp_id,
      px = point$x,
      py = point$y,
      pval = current_vals$present,
      pvalf = current_vals$future,
      envelope_min = current_vals$envelope_min,
      envelope_max = current_vals$envelope_max,
      fx = if(is.null(future_xy)) NA_real_ else future_xy[1],
      fy = if(is.null(future_xy)) NA_real_ else future_xy[2],
      fval = if(is.null(future_val)) NA_real_ else future_val,
      distance = if(is.null(distance)) NA_real_ else distance,
      directionr = if(is.null(direction)) NA_real_ else direction,
      directiondeg = if(is.null(direction)) NA_real_ else (450 - direction * (180/pi)) %% 360,
      area = area
    )
  }
  
  # Process chunk of points
  process_chunk <- function(chunk_info) {
    chunk_start <- chunk_info$start
    chunk_end <- chunk_info$end
    chunk_number <- chunk_info$chunk_num 
    
    # Initialize results list for SST and pH
    results_list <- list(
      data.table(matrix(NA_real_, nrow=0, ncol=length(result_columns))),
      data.table(matrix(NA_real_, nrow=0, ncol=length(result_columns)))
    )
    setnames(results_list[[1]], result_columns)
    setnames(results_list[[2]], result_columns)
    
    # Unwrap mask and future rasters
    mask <- unwrap(mask_wrapped)
    temp_mask <- mask
    future_temp <- unwrap(future_temp_wrapped)
    future_ph <- unwrap(future_ph_wrapped)
    
    # Process each point in the chunk
    for(point_idx in chunk_start:chunk_end) {
      current_point <- unique_points[point_idx]
      
      # Get species at this point efficiently
      species_at_point <- unique(sp_coords[x == current_point$x & 
                                             y == current_point$y, species_id])
      
      # Get all species data for this point
      point_data <- sp_coords[x == current_point$x & y == current_point$y]
      
      # Get future values for this point
      point_coords <- matrix(c(current_point$x, current_point$y), ncol=2)
      future_temp_val <- terra::extract(future_temp, point_coords)[[1]]
      future_ph_val <- terra::extract(future_ph, point_coords)[[1]]
      
      #calculate distance raster
      temp_mask[] <- mask[]
      origin_cell <- cellFromXY(temp_mask, point_coords)
      temp_mask[origin_cell] <- 500
      dist_raster <- gridDist(temp_mask, target=500, scale=1000)

      # Process each unique species at this point
      for(sp_id in species_at_point) {
        # Get envelope for current species
        species_envelope <- sp_envelope[.(sp_id)]
        
        # Get data for current species at this point
        species_point_data <- point_data[point_data$species_id == sp_id,][1]
        
        if(nrow(species_envelope) == 0) {
          # Create empty current_vals structure for both SST and pH
          sst_vals <- list(
            present = species_point_data$sst_mean,
            future = future_temp_val,
            envelope_min = NA_real_,
            envelope_max = NA_real_
          )
          ph_vals <- list(
            present = species_point_data$ph_mean,
            future = future_ph_val,
            envelope_min = NA_real_,
            envelope_max = NA_real_
          )
          # Add NA rows for both SST and pH when no envelope exists
          results_list[[1]] <- rbindlist(list(results_list[[1]], 
                                              create_result_row(current_point, sp_id, sst_vals, point_data$area[1])),fill=T)
          results_list[[2]] <- rbindlist(list(results_list[[2]], 
                                              create_result_row(current_point, sp_id, ph_vals, point_data$area[1])),fill=T)
          next
        }
        
      # Process SST and pH
      for(i in 1:2) {  # 1 = SST, 2 = pH
        var_name <- if(i == 1) "sst_future" else "ph_future"
        future_raster <- if(i == 1) future_temp else future_ph
        future_val <- if(i == 1) future_temp_val else future_ph_val
        # Set up current values and envelope bounds
        current_vals <- if(i == 1) {
          list(
            present = species_point_data$sst_mean,
            future = future_temp_val,
            envelope_min = species_envelope$sst_min,
            envelope_max = species_envelope$sst_max
          )
        } else {
          list(
            present = species_point_data$ph_mean,
            future = future_ph_val,
            envelope_min = species_envelope$ph_min,
            envelope_max = species_envelope$ph_max
          )
        }
        
        # Skip if missing values
        if(any(is.na(unlist(current_vals)))) {
          results_list[[i]] <- rbindlist(list(results_list[[i]], 
                                              create_result_row(current_point, sp_id, current_vals, point_data$area[1])), fill=TRUE)
          next
        }
        
        # Check if future value at origin point is within envelope
        if(!is.na(future_val) && 
           future_val >= current_vals$envelope_min && 
           future_val <= current_vals$envelope_max) {
          # Future value at origin is within envelope, create row with distance 0
          results_list[[i]] <- rbindlist(list(results_list[[i]], 
                                              create_result_row(current_point, sp_id, current_vals, point_data$area[1],
                                                                future_xy = c(current_point$x, current_point$y),
                                                                future_val = future_val,
                                                                distance = 0,
                                                                direction = 0)), fill=TRUE)
          next
        }
        
        # if origin point is not within envelope, find cells within envelope range
        matches <- which(!is.na(future_raster[]) & 
                           future_raster[] <= current_vals$envelope_max & 
                           future_raster[] >= current_vals$envelope_min)
          
          if(length(matches)>0) {
            # Calculate distances to matches
            matches_xy <- xyFromCell(future_raster, matches)
            distances <- extract(dist_raster, matches_xy)$landseamask
            valid_distances <- !is.na(distances)
            
            if(any(valid_distances)) {
              min_distance <- min(distances[valid_distances])
              closest_idx <- which(distances == min_distance)
              
              for(k in closest_idx) {
                # Calculate direction
                direction <- atan2(matches_xy[k,2] - current_point[["y"]],
                                   matches_xy[k,1] - current_point[["x"]])

                # Get the future conditions at the closest matched locations
                matched_val <- terra::extract(future_raster,matrix(matches_xy[k,], ncol=2))[[1]]
                
                # Create new row
                results_list[[i]] <- rbindlist(list(results_list[[i]], 
                                                    create_result_row(current_point, sp_id, current_vals, point_data$area[1],
                                                                      future_xy = matches_xy[k,],
                                                                      future_val = matched_val,
                                                                      distance = distances[k],
                                                                      direction = direction)), fill=TRUE)
                }
              } else {
                # All distances were NA
                results_list[[i]] <- rbindlist(list(results_list[[i]], 
                                                    create_result_row(current_point, sp_id, current_vals, point_data$area[1])), fill=TRUE)
            }
            } else {
              # No matches within envelope 
              results_list[[i]] <- rbindlist(list(results_list[[i]], 
                                                  create_result_row(current_point, sp_id, current_vals, point_data$area[1])), fill=TRUE)
        }
      }
      }
      # Garbage collection every 100 points
      if(point_idx %% 250  == 0) {
        rm(dist_raster)
        gc(verbose = FALSE)
      }
    }
    
    return(results_list)
  }
  
  # Process chunks
  chunks <- lapply(1:n_chunks, function(i) {
    start_idx <- (i-1) * points_per_chunk + 1
    end_idx <- min(i * points_per_chunk, total_points)
    list(start = start_idx, end = end_idx, chunk_num = i)
  })
  
  results <- parLapply(cl, chunks, process_chunk)
  
  stopCluster(cl)
  
  # Combine results
  cat("\nCombining results...\n")
  final_results <- list(
    sst = rbindlist(lapply(results, function(x) x[[1]]), fill=TRUE),
    ph = rbindlist(lapply(results, function(x) x[[2]]), fill=TRUE)
  )
  
  expected_rows <- nrow(sp_coords)
  actual_rows_sst <- nrow(final_results$sst)
  actual_rows_ph <- nrow(final_results$ph)
  cat(sprintf("Expected %d rows, got %d sst rows and %d ph rows\n", 
              expected_rows, actual_rows_sst,actual_rows_ph))
  
  end_time <- Sys.time()
  total_time <- as.numeric(end_time - start_time, units="secs")
  cat(sprintf("\nProcessing completed in %0.1f seconds (%0.1f minutes)\n", 
              total_time, total_time/60))

  return(final_results)
}




# Modified analogdifference for envelopes, run per origin location (is slower than analogdifference, but same calculation)
analogdifference_origin <- function(tempdf, phdf, mask, future_temp, future_ph, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  start_time <- Sys.time()
  
  # Convert inputs to data.table if they aren't already
  if (!is.data.table(tempdf)) tempdf <- as.data.table(tempdf)
  if (!is.data.table(phdf)) phdf <- as.data.table(phdf)
  
  # Pre-process future values into data.tables for faster lookup
  cat("ANALOGDIFF: Pre-processing future values...\n")
  future_rasters <- c(future_temp,future_ph)
  future_val_dt <- as.data.table(terra::as.data.frame(future_rasters, xy = TRUE))
  setnames(future_val_dt, c("x", "y", "future_temp","future_ph"))
  setkey(future_val_dt, x, y)

  # Create database of unique coordinates
  temp_coords <- tempdf[, .(px = px, py = py)]
  ph_coords <- phdf[, .(px = px, py = py)]
  all_coords <- unique(rbind(temp_coords, ph_coords))
  setkey(all_coords, px, py)
  all_coords[, point_id := .I]
  
  # Get unique species IDs
  unique_species <- unique(c(tempdf$species_id, phdf$species_id))
  cat(sprintf("Processing %d unique coordinates for %d species...\n", 
              nrow(all_coords), length(unique_species)))
  
  # Set number of cores and chunk size
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  total_points <- nrow(all_coords)
  chunk_size <- min(5000,max(250, floor(total_points/(4 * ncores))))
  n_chunks <- ceiling(total_points/chunk_size)
  
  # Create wrapped version of mask for parallel processing
  mask_wrapped <- wrap(mask)
  
  cat(sprintf("Processing %d unique points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("all_coords", "tempdf", "phdf", "mask_wrapped","future_val_dt"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
    setkey(future_val_dt, x, y)
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
      point_data_temp <- tempdf[px == current_point$px & py == current_point$py]
      point_data_ph <- phdf[px == current_point$px & py == current_point$py]
      
      # Process each unique species at this point
      unique_species_at_point <- unique(c(point_data_temp$species_id, point_data_ph$species_id))
      
      for(sp_id in unique_species_at_point) {
        
        # Get temperature and pH data for current location
        temp_data <- tempdf[.(species_id = sp_id, px = current_point$px, py = current_point$py), on = .(species_id, px, py)]
        
        ph_data <- phdf[.(species_id = sp_id, px = current_point$px, py = current_point$py), on = .(species_id, px, py)]
        
        # Skip if either dataset is missing this location
        if(nrow(temp_data) == 0 || nrow(ph_data) == 0 ||
           all(is.na(temp_data$fx)) || all(is.na(ph_data$fx))) {
          new_row <- create_na_row(current_point, temp_data, ph_data, sp_id)
          results <- rbindlist(list(results, new_row))
          next
        }
        
        # Check if both analogs are at origin (either by coordinates or distance)
        is_temp_at_origin <- any(!is.na(temp_data$fx) & 
                                   !is.na(temp_data$distance) &
                                   ((temp_data$fx == current_point$px & temp_data$fy == current_point$py) |
                                      temp_data$distance == 0))
        
        is_ph_at_origin <- any(!is.na(ph_data$fx) & 
                                 !is.na(ph_data$distance) &
                                 ((ph_data$fx == current_point$px & ph_data$fy == current_point$py) |
                                    ph_data$distance == 0))
        
        if(is_temp_at_origin && is_ph_at_origin) {
          # Get the rows where analogs are at origin
          temp_origin_row <- temp_data[!is.na(fx) & !is.na(distance) & 
                                         ((fx == current_point$px & fy == current_point$py) |
                                            distance == 0)][1]
          ph_origin_row <- ph_data[!is.na(fx) & !is.na(distance) & 
                                     ((fx == current_point$px & fy == current_point$py) |
                                        distance == 0)][1]
          
          results <- rbindlist(list(results, create_result_row(
            temp_origin_row, ph_origin_row, 0,
            temp_origin_row$pvalf, ph_origin_row$pvalf, sp_id
          )))
          next
        }
 
        min_distance <- Inf
        best_combo <- NULL
        
        # Calculate distances between temperature and pH analogs
        for(i in seq_len(nrow(temp_data))) {
          if(is.na(temp_data$fx[i])) {
            # Create NA row when skipping due to missing temperature analog
            if(i == nrow(temp_data)) {  # Only if this is the last iteration and no valid combo found
              if(is.null(best_combo)) {
                new_row <- create_na_row(current_point, temp_data, ph_data, sp_id)
                results <- rbindlist(list(results, new_row))
              }
            }
            next
          }
          
          # Reset temp_mask
          temp_mask[] <- mask[]
          tempf_cell <- cellFromXY(temp_mask, matrix(c(temp_data$fx[i], temp_data$fy[i]), ncol=2))
          temp_mask[tempf_cell] <- 500
          
          dist_raster <- gridDist(temp_mask, target=500, scale=1000)
          
          for(j in seq_len(nrow(ph_data))) {
            if(is.na(ph_data$fx[j])) {
              # Create NA row when skipping due to missing temperature analog
              if(j == nrow(ph_data)) {  # Only if this is the last iteration and no valid combo found
                if(is.null(best_combo)) {
                  new_row <- create_na_row(current_point, temp_data, ph_data, sp_id)
                  results <- rbindlist(list(results, new_row))
                }
              }
              next
            }
            
            distance <- extract(dist_raster,matrix(c(ph_data$fx[j], ph_data$fy[j]), ncol=2))$landseamask
            
            if(!is.na(distance) && distance <= min_distance) {
              min_distance <- distance
              
              # Extract future pH at temperature analog location
              fanalogph_temp <- future_val_dt[.(x=temp_data$fx[i], y=temp_data$fy[i]), future_ph]
              
              # Extract future temperature at pH analog location
              fanalogtemp_ph <- future_val_dt[.(x=ph_data$fx[j], y=ph_data$fy[j]), future_temp]                            
              
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
      if(point_idx %% 250 == 0) {
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
      px_temp = point$px,
      py_temp = point$py,
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
      px_ph = point$px,
      py_ph = point$py,
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

#analogdiff calculation run per ph analog
analogdifference <- function(tempdf, phdf, mask, future_temp, future_ph, ncores = NULL) {
  require(parallel)
  require(data.table)
  require(terra)
  
  start_time <- Sys.time()
  
  # Convert inputs to data.table if they aren't already
  if (!is.data.table(tempdf)) tempdf <- as.data.table(tempdf)
  if (!is.data.table(phdf)) phdf <- as.data.table(phdf)
  
  # Pre-process future values into data.tables for faster lookup
  cat("ANALOGDIFF: Pre-processing future values...\n")
  future_rasters <- c(future_temp, future_ph)
  future_val_dt <- as.data.table(terra::as.data.frame(future_rasters, xy = TRUE))
  setnames(future_val_dt, c("x", "y", "future_temp", "future_ph"))
  setkey(future_val_dt, x, y)
  
  # Create merged dataset of all possible combinations
  cat("Creating merged dataset of all combinations...\n")
  merged_dt <- merge(tempdf, phdf, 
                     by = c("species_id", "px", "py"),
                     suffixes = c("_temp", "_ph"),
                     allow.cartesian = TRUE)
  
  # Split into cases with and without analogs
  cases_with_analogs <- merged_dt[!is.na(fx_temp) & !is.na(fx_ph)]
  cases_without_analogs <- merged_dt[is.na(fx_temp) | is.na(fx_ph)]
  
  # Get unique pH analog locations
  unique_ph_analogs <- unique(cases_with_analogs[, .(fx_ph, fy_ph)])
  cat(sprintf("Processing %d unique pH analog locations...\n", nrow(unique_ph_analogs)))
  
  # Set up parallel processing
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  chunk_size <- min(1000, max(100, floor(nrow(unique_ph_analogs)/ncores)))
  n_chunks <- ceiling(nrow(unique_ph_analogs)/chunk_size)
  
  # Create wrapped version of mask for parallel processing
  mask_wrapped <- wrap(mask)
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("cases_with_analogs", "mask_wrapped", "future_val_dt"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
    setkey(future_val_dt, x, y)
  })
  
  # Helper function to create NA row
  create_na_row <- function(data) {
    data.table(
      species_id = data$species_id,
      px_temp = data$px,
      py_temp = data$py,
      pval_temp = data$pval_temp,
      pvalf_temp = data$pvalf_temp,
      pmin_temp = data$envelope_min_temp,
      pmax_temp = data$envelope_max_temp,
      fx_temp = data$fx_temp,
      fy_temp = data$fy_temp,
      fval_temp = data$fval_temp,
      distance_temp = data$distance_temp,
      directionr_temp = data$directionr_temp,
      directiondeg_temp = data$directiondeg_temp,
      px_ph = data$px,
      py_ph = data$py,
      pval_ph = data$pval_ph,
      pvalf_ph = data$pvalf_ph,
      pmin_ph = data$envelope_min_ph,
      pmax_ph = data$envelope_max_ph,
      fx_ph = data$fx_ph,
      fy_ph = data$fy_ph,
      fval_ph = data$fval_ph,
      distance_ph = data$distance_ph,
      directionr_ph = data$directionr_ph,
      directiondeg_ph = data$directiondeg_ph,
      area = data$area_temp,
      analog_diff = NA_real_,
      fanalogph_temp = NA_real_,
      fanalogphdiff_temp = NA_real_,
      fanalogtemp_ph = NA_real_,
      fanalogtempdiff_ph = NA_real_
    )
  }
  
  # Helper function to create result row
  create_result_row <- function(data, distance) {
    data.table(
      species_id = data$species_id,
      px_temp = data$px,
      py_temp = data$py,
      pval_temp = data$pval_temp,
      pvalf_temp = data$pvalf_temp,
      pmin_temp = data$envelope_min_temp,
      pmax_temp = data$envelope_max_temp,
      fx_temp = data$fx_temp,
      fy_temp = data$fy_temp,
      fval_temp = data$fval_temp,
      distance_temp = data$distance_temp,
      directionr_temp = data$directionr_temp,
      directiondeg_temp = data$directiondeg_temp,
      px_ph = data$px,
      py_ph = data$py,
      pval_ph = data$pval_ph,
      pvalf_ph = data$pvalf_ph,
      pmin_ph = data$envelope_min_ph,
      pmax_ph = data$envelope_max_ph,
      fx_ph = data$fx_ph,
      fy_ph = data$fy_ph,
      fval_ph = data$fval_ph,
      distance_ph = data$distance_ph,
      directionr_ph = data$directionr_ph,
      directiondeg_ph = data$directiondeg_ph,
      area = data$area_temp,
      analog_diff = data$analog_diff,
      fanalogph_temp = data$fanalogph_temp,
      fanalogphdiff_temp = data$fanalogphdiff_temp,
      fanalogtemp_ph = data$fanalogtemp_ph,
      fanalogtempdiff_ph = data$fanalogtempdiff_ph
    )
  }
  
  # Process chunk of pH analogs
  process_chunk <- function(chunk_info) {
    chunk_start <- chunk_info$start
    chunk_end <- chunk_info$end
    
    # Initialize results
    results <- data.table()
    
    # Unwrap mask once at the start
    mask <- unwrap(mask_wrapped)
    
    # Process each pH analog in chunk
    for(i in chunk_start:chunk_end) {
      ph_analog <- unique_ph_analogs[i]
      
      # Get all cases that use this pH analog
      ph_cases <- cases_with_analogs[fx_ph == ph_analog$fx_ph & 
                                       fy_ph == ph_analog$fy_ph]
      
      # Skip if no cases use this analog
      if(nrow(ph_cases) == 0) next
      
      # Create distance raster for this pH analog
      temp_mask <- mask
      ph_cell <- cellFromXY(temp_mask, 
                            matrix(c(ph_analog$fx_ph, ph_analog$fy_ph), ncol=2))
      temp_mask[ph_cell] <- 500
      dist_raster <- gridDist(temp_mask, target=500, scale=1000)
      
      # Calculate distances for all temperature analogs
      ph_cases[, analog_diff := extract(dist_raster, 
                                        matrix(c(fx_temp, fy_temp), ncol=2))$landseamask]
      
      # Extract future values for each combination
      ph_cases[, `:=`(
        fanalogph_temp = future_val_dt[.(x=fx_temp, y=fy_temp), future_ph], #future ph at temp analog
        fanalogtemp_ph = future_val_dt[.(x=fx_ph, y=fy_ph), future_temp], #future temp at ph analog
        fanalogphdiff_temp = future_val_dt[.(x=fx_temp, y=fy_temp), future_ph] - pval_ph, #change in ph at temp analog
        fanalogtempdiff_ph = future_val_dt[.(x=fx_ph, y=fy_ph), future_temp] - pval_temp #change in temp at ph analog
      )]
      
      # Add to results
      formatted_cases <- ph_cases[, create_result_row(.SD, analog_diff), by = 1:nrow(ph_cases)]
      results <- rbindlist(list(results, formatted_cases))
      
      # Garbage collection
      if(i %% 100 == 0) {
        rm(dist_raster)
        gc(verbose = FALSE)
      }
    }
    
    return(results)
  }
  
  # Process chunks in parallel
  chunks <- lapply(1:n_chunks, function(i) {
    list(start = (i-1) * chunk_size + 1,
         end = min(i * chunk_size, nrow(unique_ph_analogs)))
  })
  
  # Get results for cases with analogs
  results_with_analogs <- rbindlist(parLapply(cl, chunks, process_chunk))
  stopCluster(cl)
  
  # Find minimum distance combinations
  cat("\nFinding minimum distance combinations...\n")
  best_combinations <- results_with_analogs[, .SD[which.min(analog_diff)], 
                                            by = .(species_id, px_temp, py_temp)]
  
  # Process cases without analogs (both temp and pH)
  cat("Processing cases without analogs...\n")
  na_results <- cases_without_analogs[, create_na_row(.SD), by = .(species_id, px_temp = px, py_temp = py, px_ph=px, py_ph=py)]
  
  # Combine results
  cat("\nCombining results...\n")
  final_results <- rbindlist(list(best_combinations, na_results),fill=T)
  
  # Sort results
  setorder(final_results, species_id, px_temp, py_temp)
  
  # Verify results
  sp_coords <- unique(paste(tempdf$species_id, tempdf$px, tempdf$py))
  expected_rows <- length(sp_coords)
  actual_rows <- nrow(final_results)
  cat(sprintf("Expected %d rows, got %d rows\n", 
              expected_rows, actual_rows))
  
  total_time <- as.numeric(Sys.time() - start_time, units="secs")
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
  
  terra::terraOptions(memfrac=0.8,progress=10)
  
  # Validate inputs
  if (!all(sapply(list(future_temp, future_ph), 
                  function(x) inherits(x, "SpatRaster")))) {
    stop("future_temp and future_ph must be SpatRaster objects")
  }
  
  # Convert sp_coords to data.table if it isn't already
  sp_coords <- as.data.table(sp_coords)
  sp_envelope <- as.data.table(sp_envelope)
  
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
  
  # find unique points & species
  unique_points <- unique(sp_coords[, .(x, y)])
  setkey(unique_points, x, y)  
  unique_points[, point_id := .I]
  unique_species <- unique(sp_coords$species_id)
  
  # Set number of cores
  if(is.null(ncores)) ncores <- max(1, parallel::detectCores() - 1)
  
  # Create wrapped version of mask for parallel processing
  cat("Creating memory-mapped mask...\n")
  mask_wrapped <- wrap(mask)
  tos_wrapped <- wrap(future_temp)
  ph_wrapped <- wrap(future_ph)
  
  # Calculate optimal chunk size
  total_points <- nrow(unique_points)
  chunk_size <- min(500, max(250, floor(total_points/20)))
  n_chunks <- ceiling(total_points/chunk_size)
  
  cat(sprintf("Processing %d unique points in %d chunks using %d cores\n", 
              total_points, n_chunks, ncores))
  
  # Initialize cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, c("unique_points", "all_data", "mask_wrapped","tos_wrapped","ph_wrapped"), 
                envir = environment())
  clusterEvalQ(cl, {
    library(terra)
    library(data.table)
  })
  
  # Helper function to create NA row
  create_result_row <- function(data, future_xy = NULL, future_vals = NULL, distance = NULL, direction = NULL) {
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
      fx = if(is.null(future_xy)) NA_real_ else future_xy[1],
      fy = if(is.null(future_xy)) NA_real_ else future_xy[2],
      fval_temp = if(is.null(future_vals)) NA_real_ else future_vals[1],
      fval_ph = if(is.null(future_vals)) NA_real_ else future_vals[2],
      distance = if(is.null(distance)) NA_real_ else distance,
      directionr = if(is.null(direction)) NA_real_ else direction,
      directiondeg = if(is.null(direction)) NA_real_ else (450 - direction * (180/pi)) %% 360
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
      current_point <- unique_points[point_idx]
      
      # Get species at this point 
      species_at_point <- unique(all_data[x == current_point$x & 
                                             y == current_point$y, species_id])
      
      # Get data for current location
      current_data <- all_data[x == current_point$x & y == current_point$y]
      
      if(nrow(current_data) == 0) {
        results <- rbind(results, create_result_row(current_data))
        next
      }
      
      #process each unique species at this point
      for(sp_id in species_at_point){
        #get data for current species
        species_data <- current_data[current_data$species_id == sp_id,][1]
        
        # Check if current location is still within both envelopes
      temp_in_envelope <- !is.na(species_data$sst_future) && 
        species_data$sst_future >= species_data$temp_min && 
        species_data$sst_future <= species_data$temp_max
      
      ph_in_envelope <- !is.na(species_data$ph_future) && 
        species_data$ph_future >= species_data$ph_min && 
        species_data$ph_future <= species_data$ph_max
      
      # If both conditions are met at origin, create result with distance 0
      if(temp_in_envelope && ph_in_envelope) {
        results <- rbind(results, create_result_row(
          data = species_data,
          future_xy = c(species_data$x, species_data$y),
          future_vals = c(species_data$sst_future, species_data$ph_future),
          distance = 0,
          direction = 0
        ))
        next
      }
      
      # if origin is not analog, create raster masks for valid regions
      temp_valid <- future_temp <= species_data$temp_max & future_temp >= species_data$temp_min
      ph_valid <- future_ph <= species_data$ph_max & future_ph >= species_data$ph_min
      
      # Combine conditions
      valid_cells <- which(!is.na(temp_valid[]) & temp_valid[] & 
                             !is.na(ph_valid[]) & ph_valid[] & 
                             !is.na(mask[]))
      
      if(length(valid_cells) == 0) {
        results <- rbind(results, create_result_row(species_data))
        next
      }
      
      # Get coordinates of valid cells
      valid_xy <- xyFromCell(future_temp, valid_cells)
      
      # Calculate distances for valid points
      temp_mask <- mask
      origin_cell <- cellFromXY(temp_mask, 
                                matrix(c(species_data$x, species_data$y), ncol=2))
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
          closest_match$y - species_data$y,
          closest_match$x - species_data$x
        )
        
        closest_temp <- terra::extract(future_temp, 
                                       matrix(c(closest_match$x, closest_match$y), ncol=2))[[1]]
        closest_ph <- terra::extract(future_ph, 
                                     matrix(c(closest_match$x, closest_match$y), ncol=2))[[1]]
        
        # Add result
        results <- rbind(results, create_result_row(
          data = species_data,
          future_xy = c(closest_match$x, closest_match$y),
          future_vals = c(closest_temp, closest_ph),
          distance = closest_match$distance,
          direction = direction
        ))
      } else {
        results <- rbind(results, create_result_row(species_data))
      }
      
      
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
  expected_rows <- nrow(all_data)
  actual_rows <- nrow(final_results)
  cat(sprintf("Expected %d rows, got %d rows\n", 
              expected_rows, actual_rows))
  
  total_time <- as.numeric(end_time - start_time, units="secs")
  cat(sprintf("\nProcessing completed in %0.1f seconds (%0.1f minutes)\n", 
              total_time, total_time/60))
  
  return(final_results)
}




# runs analogcalc, analogdifference, dualanalog calculations
run_analog_analysis <- function(species_list, species_envelopes, 
                                future_temp, future_ph, 
                                time_period, base_period="hist", mask, output_dir,
                                scenario, model, ncores=NULL) {
  require(terra)
  require(data.table)
  
  # Validate base_period
  if(!base_period %in% c("hist", "1985")) {
    stop("base_period must be either 'hist' or '1985'")
  }
  
  # Convert inputs to data.tables if not already
  species_dt <- as.data.table(species_list)
  envelopes_dt <- as.data.table(species_envelopes)
  
  cat("Processing all coordinates...\n")
  
  # Create output filenames
  output_files <- list(
    dual = file.path(output_dir, 
                     sprintf("dualanalog_%s_%s_%s-%s.csv",
                             model, scenario, base_period, time_period)),
    tos = file.path(output_dir,
                    sprintf("analog_tos_%s_%s_%s-%s.csv",
                            model, scenario, base_period, time_period)),
    ph = file.path(output_dir,
                   sprintf("analog_ph_%s_%s_%s-%s.csv",
                           model, scenario, base_period, time_period)),
    diff = file.path(output_dir,
                     sprintf("analogdiff_%s_%s_%s-%s.csv",
                             model, scenario, base_period, time_period))
  )
  
  # Check if output files exist
  all_files_exist <- all(file.exists(c(output_files$dual, output_files$tos, output_files$ph, output_files$diff)))
  
  if(all_files_exist) {
    cat(sprintf("\nAll output files already exist for %s to %d. Skipping analysis...\n", 
                time_period, end_year))
    return(invisible(NULL))
  }
  
  # Run dual analog analysis
  if(!file.exists(output_files$dual)) {
    cat(sprintf("\nRunning dual analog analysis for %s to %d...\n", base_period, time_period))
    dual_results <- dualanalog(
      sp_coords = species_dt,
      sp_envelope = envelopes_dt,
      future_temp = future_temp,
      future_ph = future_ph,
      mask = mask
    )
    cat("\nSaving dual analog results...\n")
    fwrite(dual_results, output_files$dual)
  } else {
    cat("\nDual analog results already exist. Skipping...\n")
  }
  
  # Run analog calc which returns a list of sst and ph results
  tosph_needed <- !file.exists(output_files$tos) || !file.exists(output_files$ph)
  if(tosph_needed) {
    cat(sprintf("\nRunning temperature and pH analog analysis for %s to %d...\n", 
                base_period, time_period))
    tosph_results <- analogcalc(
      sp_coords = species_dt,
      future_temp = future_temp,
      future_ph = future_ph,
      sp_envelope = envelopes_dt,
      mask = mask
    )
  
    # Save individual results if they don't exist
    if(!file.exists(output_files$tos)) {
      cat("\nSaving temperature analog results...\n")
      fwrite(tosph_results$sst, output_files$tos, row.names = FALSE)
    }
    if(!file.exists(output_files$ph)) {
      cat("\nSaving pH analog results...\n")
      fwrite(tosph_results$ph, output_files$ph, row.names = FALSE)
    }
  } else {
    cat("\nTemperature and pH analog results already exist. Skipping...\n")
  }
  
  # Run analog difference
  if(!file.exists(output_files$diff)) {
    cat("\nRunning analog difference analysis...\n")
    # If we didn't run tosph analysis, we need to read the files
    if(!tosph_needed) {
      tosph_results <- list(
        sst = fread(output_files$tos),
        ph = fread(output_files$ph)
      )
    }
    
    diff_results <- analogdifference(
      tempdf = tosph_results$sst,
      phdf = tosph_results$ph,
      mask = mask,
      future_temp = future_temp,
      future_ph = future_ph,
      ncores = ncores
    )
    
    cat("\nSaving analog difference results...\n")
    fwrite(diff_results, output_files$diff, row.names = FALSE)
  } else {
    cat("\nAnalog difference results already exist. Skipping...\n")
  }
  
  # Clean up
  gc()
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




# threshold functions

#function to change mon timesteps to summer-only months(summerN, summerS) then to annual summer mean 
#summerN = July-September (N Hemisphere), summerS = January-March (S Hemisphere)
#note: requires monthly input, no yr input. modelmedian input
#variable = 'tos' or 'ph'
#scenario = "hist","ssp126","ssp245","ssp370","ssp585"
#model = 'modelmean' or 'modelmedian
summerannual <- function(variable, scenario, model) {
  files <- list.files(paste0(cmip_folder,variable), 
                      pattern = paste0(variable,".+Omon.",model,".+",scenario),
                      full.names = T)
  
  timestepin <- substring(sub(".*_O","",files[1]),1,2)
  raster <- rast(files)
  
  fileoutN <- paste(sub("_Omon","_OyrsummerN",files[1]))
  fileoutS <- paste(sub("_Omon","_OyrsummerS",files[1]))
  
  # for ssp126,ssp245,ssp370,ssp585
  # 1032 months total, 86 years
  years <- 86
  
  # northern hemisphere summer (July-September)
  julyaugsept <- c(7,8,9)
  julyaugsepts <- numeric(0)
  for (i in 0:(years - 1)) {
    julyaugsepts <- c(julyaugsepts, julyaugsept + i * 12)
  }
  summerN <- subset(raster, julyaugsepts)
  
  # southern hemisphere summer (January-March)
  janfebmarch <- c(1,2,3)
  janfebmarchs <- numeric(0)
  for (i in 0:(years - 1)) {
    janfebmarchs <- c(janfebmarchs, janfebmarch + i * 12)
  }
  summerS <- subset(raster, janfebmarchs)
  
  yearIndex <- rep(1:86, each=3)  # 3 months per summer, for 86 years

  # Calculate annual summer means
  annualN <- tapp(summerN, index=yearIndex, fun=mean, na.rm=T, cores=18, filename=fileoutN)
  
  annualS <- tapp(summerS, index=yearIndex, fun=mean, na.rm=T, cores=18, filename=fileoutS)
  
  annualNS <- list(annualN, annualS)
  
  return(annualNS)
}



# dhw functions

# calculate monthly DHW 
calculate_dhw <- function(sst_file, mmm){
  sst <- rast(sst_file)
  anomalies <- sst - mmm
  anomalies[anomalies < 1] <- 0
  
  n_layers <- nlyr(anomalies)
  dhw_list <- vector("list", n_layers)
  
  for(i in 1:n_layers) {
    # get 3 month window (equivalent to 12 weeks)
    window_size <- min(3, i)
    start_idx <- max(1, i - window_size + 1)
    window <- anomalies[[start_idx:i]]
    
    # multiply by 4.34 to convert monthly anomalies to DHW
    dhw_list[[i]] <- sum(window, na.rm=TRUE) * 4.34
  }
  
  dhw_stack <- rast(dhw_list)
  return(dhw_stack)
}

# Calculate annual maximum DHW (capturing summer bleaching events)
calculate_annual_max_dhw <- function(dhw_stack){
  n_years <- nlyr(dhw_stack) %/% 12
  annual_max_dhw <- rast()
  
  for(year in 1:n_years) {
    year_start <- (year-1)*12 + 1
    year_end <- year*12
    year_data <- dhw_stack[[year_start:year_end]]
    max_dhw <- max(year_data, na.rm=T)
    annual_max_dhw <- c(annual_max_dhw, max_dhw)
  }
  
  return(annual_max_dhw)
}

# year annual severe bleaching begins 
# (after this year, there are severe bleaching events of DHW >8 every subsequent year)
year_ASB <- function(annual_max_dhw){
  n_years <- nlyr(annual_max_dhw)
  years <- 2015:2100
  
  asb_yr <- app(annual_max_dhw, function(x) {
    if(all(is.na(x))) return(NA)
    
    for(i in 1:length(x)) {
      if(all(x[i:length(x)] > 8, na.rm=TRUE)) {
        return(years[i])
      }
    }
    return(NA)
  })
  
  return(asb_yr)
}

# length of bleaching events per year
# (how many months per year are potential bleaching events, where DHW > 8)
asb_duration_annual <- function(dhw_stack){
  # mask for severe bleaching events
  severe_bleaching <- dhw_stack > 8
  
  n_years <- 86
  year_index <- rep(1:n_years, each=12)
  
  asb_duration <- tapp(severe_bleaching, 
                       index = year_index, 
                       fun = sum, 
                       na.rm = TRUE)
  
  return(asb_duration)
}


# Calculate decadal mean 
# use on output of calculate_annual_max_dhw() to find decadal max DHW 
# use on output of asb_duration_annual() to find decadal asb duration
calc_decadal_from_annual <- function(annual_dhw_raster){
  
  decade_index <- rep(1:9, each=10)[1:nlyr(annual_dhw_raster)]
  
  # Calculate decadal means
  decadal_means <- tapp(annual_dhw_raster, 
                        index=decade_index,
                        fun=mean,
                        na.rm=TRUE)
  
  return(decadal_means)
}


