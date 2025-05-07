


# Aragonite threshold exceedance analysis
# Using decadal data files (9 time points from 2020 to 2100)

scenarios <- c("ssp126","ssp245","ssp370","ssp585")
decades <- c(2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100)

# Convert scenario format for aragonite files
scenario_map <- c("126", "245", "370", "585")
names(scenario_map) <- scenarios

for(scenario in scenarios) {
  scen_file <- scenario_map[scenario]
  
  # Load aragonite data for this scenario
  # Note: this data is decadal rather than annual
  arag_file <- paste0(cmip_folder, "arag/Aragonite_median_", scen_file, "_foc.nc")
  arag_raster <- rast(arag_file)
  names(arag_raster) <- as.character(decades)
  
  # Load species data with aragonite information
  species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992_with_arag.csv'))
  species_coord <- fread(paste0(output_directory, 'species_pixels_1982-1992_with_arag.csv'))
  species_coord$x <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
  setkey(species_coord, x, y)
  
  # Get unique coordinates
  coords_dt <- unique(species_coord[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))
  
  # Calculate pixel areas
  area_raster <- cellSize(arag_raster[[1]], unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  # Merge with species data
  species_coord <- merge(species_coord, coords_dt, by = c("x", "y"))
  species_coord[, sp_area := pixel_area * pixel_cover]
  setkey(species_coord, species_id)
  
  # Extract aragonite values for each point
  arag_values <- as.data.table(terra::extract(arag_raster, points, ID=F))
  names(arag_values) <- as.character(decades)
  
  # Combine coordinates with aragonite values
  climate_dt <- copy(coords_dt)
  for(i in 1:length(decades)) {
    decade <- decades[i]
    climate_dt[, paste0("arag_", decade) := arag_values[[i]]]
  }
  
  all_species_summary <- data.table()
  
  for(sp in unique(species_envelope$id_no)) {
    sp_envelope <- species_envelope[species_envelope$id_no == sp,]
    
    sp_summary <- data.table(
      species_id = sp,
      scenario = scenario,
      decade = decades,
      mean_arag_exceed = numeric(length(decades)),     # mean of only exceeding values
      mean_arag_diff = numeric(length(decades)),       # mean difference from envelope
      arag_exceed_area = numeric(length(decades)),
      total_area = sum(species_coord[species_id == sp, sp_area])
    )
    
    # Get coords for just this species
    sp_coords <- species_coord[species_id == sp]
    
    sp_pixels_all_decades <- data.table()
    
    for(decade in decades) {
      decade_chr <- as.character(decade)
      
      # Get current decade's aragonite values
      decade_climate <- climate_dt[, .(
        x = x,
        y = y,
        arag = get(paste0("arag_", decade_chr))
      )]
      
      # Get species coordinates and current decade's aragonite values
      sp_pixels <- merge(sp_coords, decade_climate, by = c("x", "y"))
      
      # Calculate exceedances compared to envelope
      sp_pixels[, `:=`(
        decade = as.numeric(decade_chr),
        scenario = scenario,
        arag_diff = arag - sp_envelope$aragonite_min,          # all differences
        arag_exceed = pmin(0, arag - sp_envelope$aragonite_min)  # only below threshold (negative values)
      )]
      
      sp_pixels_all_decades <- rbindlist(list(sp_pixels_all_decades, 
                                              sp_pixels[, .(x, y, decade, scenario, 
                                                            species_id = sp,
                                                            arag_diff,
                                                            arag_exceed,
                                                            pixel_area, 
                                                            pixel_cover,
                                                            sp_area)]))
      
      # Calculate summary stats for this decade
      decade_idx <- which(decades == decade)
      sp_summary[decade_idx, `:=`(
        # Mean exceedance (only for exceeding pixels)
        mean_arag_exceed = weighted.mean(sp_pixels[arag_exceed < 0, arag_exceed], 
                                         w = sp_pixels[arag_exceed < 0, sp_area], 
                                         na.rm = TRUE),
        # Mean difference from envelope (all pixels)
        mean_arag_diff = weighted.mean(sp_pixels$arag_diff, 
                                       w = sp_pixels$sp_area, 
                                       na.rm = TRUE),
        # Areas of exceedance
        arag_exceed_area = sum(sp_pixels[arag_exceed < 0, sp_area], 
                               na.rm = TRUE),
        total_area = sum(sp_pixels$sp_area, 
                         na.rm = TRUE)
      )]
    }
    
    # Save pixel-level results for this species
    fwrite(sp_pixels_all_decades, 
           file = paste0(output_directory, "threshold/species_", sp, "_exceedance_arag_", scenario, ".csv"))
    
    # Add to summary
    all_species_summary <- rbindlist(list(all_species_summary, sp_summary), fill=TRUE)
  }
  
  # Save summary file
  fwrite(all_species_summary, 
         file = paste0(output_directory, "species_exceedance_summary_arag_", scenario, ".csv"))
}












# Run analysis for all coral
for(scenario in scenarios) {
  scen_file <- scenario_map[scenario]
  
  # Load aragonite data
  arag_file <- paste0(cmip_folder, "arag/Aragonite_median_", scen_file, "_foc.nc")
  arag_raster <- rast(arag_file)
  names(arag_raster) <- as.character(decades)
  
  # Load allcoral data
  allcoral_envelope <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992_with_arag.csv'))
  allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992_with_arag.csv'))
  allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)
  
  # Get unique coordinates
  coords_dt <- unique(allcoral_coord[, .(x, y, pixel_cover)]) 
  setkey(coords_dt, x, y)
  points <- vect(coords_dt[, .(x, y)], geom = c("x", "y"))  
  
  # Calculate pixel areas
  area_raster <- cellSize(arag_raster[[1]], unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  # Calculate species area based on pixel coverage
  coords_dt[, sp_area := pixel_area * pixel_cover]  
  
  # Extract aragonite values
  arag_values <- as.data.table(terra::extract(arag_raster, points, ID=F))
  names(arag_values) <- as.character(decades)
  
  # Combine coordinates with aragonite values
  climate_dt <- copy(coords_dt)
  for(i in 1:length(decades)) {
    decade <- decades[i]
    climate_dt[, paste0("arag_", decade) := arag_values[[i]]]
  }
  
  allcoral_summary <- data.table(
    species_id = "allcoral",
    scenario = scenario,
    decade = decades,
    mean_arag_exceed = numeric(length(decades)),
    mean_arag_diff = numeric(length(decades)),
    arag_exceed_area = numeric(length(decades)),
    total_area = sum(coords_dt$sp_area, na.rm = TRUE)  
  )
  
  allcoral_pixels_all_decades <- data.table()
  
  for(decade in decades) {
    decade_chr <- as.character(decade)
    
    # Get current decade's aragonite values
    decade_climate <- climate_dt[, .(
      x = x,
      y = y,
      pixel_area = pixel_area,
      pixel_cover = pixel_cover,
      sp_area = sp_area,  
      arag = get(paste0("arag_", decade_chr))
    )]
    
    # Merge coordinates with aragonite data - only merge by x,y since other values are already included
    allcoral_pixels <- merge(allcoral_coord[, .(x, y, species_id="allcoral")], decade_climate, by = c("x", "y"))
    
    # Calculate exceedances
    allcoral_pixels[, `:=`(
      decade = as.numeric(decade_chr),
      scenario = scenario,
      arag_diff = arag - allcoral_envelope$aragonite_min,
      arag_exceed = pmin(0, arag - allcoral_envelope$aragonite_min)
    )]
    
    # Store pixel-level data
    allcoral_pixels_all_decades <- rbindlist(list(allcoral_pixels_all_decades, 
                                                  allcoral_pixels[, .(x, y, decade, scenario, 
                                                                      species_id,
                                                                      arag_diff,
                                                                      arag_exceed,
                                                                      pixel_area, 
                                                                      pixel_cover,
                                                                      sp_area)]))
    
    # Calculate summary statistics with proper error handling
    decade_idx <- which(decades == decade)
    
    # For mean calculations, ensure we have data and handle NAs properly
    mean_exceed <- tryCatch({
      if(sum(allcoral_pixels$arag_exceed < 0, na.rm = TRUE) > 0) {
        weighted.mean(allcoral_pixels[arag_exceed < 0, arag_exceed], 
                      w = allcoral_pixels[arag_exceed < 0, sp_area], 
                      na.rm = TRUE)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    mean_diff <- tryCatch({
      weighted.mean(allcoral_pixels$arag_diff, 
                    w = allcoral_pixels$sp_area, 
                    na.rm = TRUE)
    }, error = function(e) NA_real_)
    
    allcoral_summary[decade_idx, `:=`(
      mean_arag_exceed = mean_exceed,
      mean_arag_diff = mean_diff,
      arag_exceed_area = sum(allcoral_pixels[arag_exceed < 0, sp_area], 
                             na.rm = TRUE),
      total_area = sum(allcoral_pixels$sp_area, 
                       na.rm = TRUE)
    )]
  }
  
  # Save all coral results
  fwrite(allcoral_pixels_all_decades, 
         file = paste0(output_directory, "threshold/allcoral_exceedance_arag_", scenario, ".csv"))
  
  # Save summary file
  fwrite(allcoral_summary, 
         file = paste0(output_directory, "allcoral_exceedance_summary_arag_", scenario, ".csv"))
}
