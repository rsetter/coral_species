# calculate projected dhw

source("functions.R")

# open mmm (maximum monthly mean)
mmm <- rast(paste0(iucn_sp_folder,'sst_foc_1982-1992_max.tif'))
mmm <- rotate(mmm,left=F)
mmm <- crop(mmm,ext(-0.5, 359.5, -90, 90))

# Process all scenarios
scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

for (scenario in scenarios) {
  # Construct filename
  sst_file <- paste0(cmip_folder,"tos/tos_Omon_modelmeanfoc_", scenario, "_201501-210012.tif")
  
  # calculate monthly DHW
  dhw <- calculate_dhw(sst_file, mmm)
  writeRaster(dhw, filename = paste0(cmip_folder,"tos/dhw_Omon_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
  # calculate annual maximum DHW
  annual_max_dhw <- calculate_annual_max_dhw(dhw)
  writeRaster(annual_max_dhw,filename = paste0(cmip_folder,"tos/dhw_Oyr_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
  # calculate decadal means
  decadal_max_dhw <- calc_decadal_from_annual(annual_max_dhw)
  writeRaster(decadal_max_dhw,filename = paste0(cmip_folder,"tos/dhw_Odec_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
  # calculate year of annual severe bleaching
  asb_year <- year_ASB(annual_max_dhw)
  writeRaster(asb_year,filename = paste0(cmip_folder,"tos/dhw_ASByear_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
  #calculate duration of annual bleaching (how many months per year are potential bleaching events, where DHW > 8)
  bleaching_duration <- asb_duration_annual(dhw)
  writeRaster(bleaching_duration,filename = paste0(cmip_folder,"tos/dhw_duration_Oyr_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
  #calculate decadal mean of annual bleaching duration
  decadal_duration <- calc_decadal_from_annual(bleaching_duration)
  writeRaster(decadal_duration,filename = paste0(cmip_folder,"tos/dhw_duration_Odec_modelmeanfoc_", scenario, "_201501-210012.tif"),overwrite=T)
  
}




### SUMMARIZE RESULTS 

for(scenario in scenarios){
 
  pixel_coverage <- fread(file.path(output_directory, 'species_pixelcover.csv'))
  pixel_coverage$x <- ifelse(pixel_coverage$x < 0, pixel_coverage$x + 360, pixel_coverage$x)
  names(pixel_coverage) <- c("x", "y", "pixelcover_pct","species_id")
  setkey(pixel_coverage, x, y)
  
  coords_dt <- unique(pixel_coverage[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))
  
  # year of asb
  asb_year <- rast(paste0(cmip_folder,"tos/dhw_ASByear_modelmeanfoc_", scenario, "_201501-210012.tif"))
  asb_values <- as.data.table(terra::extract(asb_year, points))
  coords_dt[, asb_year := asb_values[[2]]]
  
  # decadal max DHW
  decadal_max_dhw <- rast(paste0(cmip_folder,"tos/dhw_Odec_modelmeanfoc_", scenario, "_201501-210012.tif"))
  dec_max_values <- as.data.table(terra::extract(decadal_max_dhw, points,ID=F))
  decades <- seq(2020, 2100, by=10)
  setnames(dec_max_values, paste0("decadal_max_dhw_", decades))
  coords_dt[, paste0("decadal_max_dhw_", decades) := dec_max_values]
  
  # decadal duration of asb
  decadal_duration <- rast(paste0(cmip_folder,"tos/dhw_duration_Odec_modelmeanfoc_", scenario, "_201501-210012.tif"))
  dec_dur_values <- as.data.table(terra::extract(decadal_duration, points,ID=F))
  setnames(dec_dur_values, paste0("decadal_duration_", decades))
  coords_dt[, paste0("decadal_duration_", decades) := dec_dur_values]
  
  # area of pixels
  area_raster <- cellSize(asb_year,unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values[[2]]]
  
  # merge with species data
  all_species_data <- merge(pixel_coverage, coords_dt, by = c("x", "y"))
  all_species_data[, sp_area := pixel_area * pixelcover_pct]
  setkey(all_species_data, species_id)
  
  
  ### summarize per species 
  
  dhw_summary_data <- all_species_data %>%
    group_by(species_id) %>%
    summarise(
      total_area = sum(sp_area, na.rm = TRUE),
      
      # ASB year statistics
      mean_asb_year = mean(asb_year, na.rm = TRUE),
      min_asb_year = min(asb_year, na.rm = TRUE),
      max_asb_year = max(asb_year,na.rm=T),
      
      # Decadal statistics
      across(starts_with("decadal_max_dhw_"),
             list(
               mean = ~weighted.mean(.x, w = sp_area, na.rm = TRUE),
               max = ~max(.x, na.rm = TRUE)
             )),
      
      across(starts_with("decadal_duration_"),
             list(
               mean = ~weighted.mean(.x, w = sp_area, na.rm = TRUE),
               max = ~max(.x, na.rm = TRUE)
             ))
    )
  
  #save summary csv
  fwrite(all_species_data,file = file.path(output_directory,paste0("coralsp_dhwsummary_detailed_", scenario, ".csv")))
  
  fwrite(dhw_summary_data, file = file.path(output_directory,paste0("coralsp_dhwsummary_", scenario, ".csv")))  
  
}















# annual DHW analysis
for(scenario in scenarios){
  # Load DHW data
  dhw_max <- rast(paste0(cmip_folder,"tos/dhw_Oyr_modelmeanfoc_", scenario, "_201501-210012.tif"))
  dhw_duration <- rast(paste0(cmip_folder,"tos/dhw_duration_Oyr_modelmeanfoc_", scenario, "_201501-210012.tif"))
  
  # Load species range data
  species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
  species_coord <- fread(file.path(output_directory, 'species_pixels_1982-1992.csv'))
  species_coord$x <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
  setkey(species_coord, x, y)
  
  # Get unique coordinates
  coords_dt <- unique(species_coord[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))
  
  # Calculate pixel areas
  area_raster <- cellSize(dhw_max, unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  # Merge with species data
  species_coord <- merge(species_coord, coords_dt, by = c("x", "y"))
  species_coord[, sp_area := pixel_area * pixel_cover]
  setkey(species_coord, species_id)
  
  # Extract values for each point
  dhw_max_values <- as.data.table(terra::extract(dhw_max, points, ID=F))
  dhw_duration_values <- as.data.table(terra::extract(dhw_duration, points, ID=F))
  
  years <- 2015:2100
  names(dhw_max_values) <- as.character(years)
  names(dhw_duration_values) <- as.character(years)
  
  climate_dt <- copy(coords_dt) 
  for(i in 1:length(years)) {
    y <- years[i]
    climate_dt[, paste0("dhw_max_", y) := dhw_max_values[[i]]]
    climate_dt[, paste0("dhw_duration_", y) := dhw_duration_values[[i]]]
  }
  
  all_species_summary <- data.table()
  
  for(sp in unique(species_envelope$id_no)){
    sp_coords <- species_coord[species_id == sp]
    
    sp_summary <- data.table(
      species_id = sp,
      scenario = scenario,
      year = years,
      mean_dhw_max = numeric(length(years)),
      mean_dhw_duration = numeric(length(years)),
      dhw_exceed_area = numeric(length(years)),      # area where DHW > 8
      duration_exceed_area = numeric(length(years)),  # area where duration > 3 months
      total_area = sum(sp_coords$sp_area)
    )
    
    for(y in years){
      y_chr <- as.character(y)
      
      # Get current year's climate values
      year_climate <- climate_dt[, .(
        x = x,
        y = y,
        dhw_max = get(paste0("dhw_max_", y_chr)),
        dhw_duration = get(paste0("dhw_duration_", y_chr))
      )]
      
      # Merge coordinates with climate data
      sp_pixels <- merge(sp_coords, year_climate, by = c("x", "y"))
      
      # Calculate exceedances
      sp_pixels[, `:=`(
        dhw_exceed = dhw_max > 8,
        duration_exceed = dhw_duration > 3
      )]
      
      # Calculate summary stats for this year
      year_idx <- which(years == y)
      sp_summary[year_idx, `:=`(
        mean_dhw_max = weighted.mean(sp_pixels$dhw_max, 
                                     w = sp_pixels$sp_area, 
                                     na.rm = TRUE),
        mean_dhw_duration = weighted.mean(sp_pixels$dhw_duration, 
                                          w = sp_pixels$sp_area, 
                                          na.rm = TRUE),
        dhw_exceed_area = sum(sp_pixels[dhw_exceed == TRUE, sp_area], 
                              na.rm = TRUE),
        duration_exceed_area = sum(sp_pixels[duration_exceed == TRUE, sp_area], 
                                   na.rm = TRUE)
      )]
    }
    
    # Add to summary
    all_species_summary <- rbindlist(list(all_species_summary, sp_summary), fill=T)
  }
  
  # Save summary file
  fwrite(all_species_summary, 
         file = paste0(output_directory, "coralsp_dhwsummary_annual_", scenario, ".csv"))
}

  # Now process all coral data
  allcoral_envelope <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
  allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))
  allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)
  
  
  # Get unique coordinates
  coords_dt <- unique(allcoral_coord[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))
  
  # Calculate pixel areas
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  allcoral_coord <- merge(allcoral_coord, coords_dt, by = c("x", "y"))
  allcoral_coord[, sp_area := pixel_area * pixel_cover]
  
  allcoral_summary <- data.table(
    species_id = "allcoral",
    scenario = scenario,
    year = years,
    mean_dhw_max = numeric(length(years)),
    mean_dhw_duration = numeric(length(years)),
    dhw_exceed_area = numeric(length(years)),
    duration_exceed_area = numeric(length(years)),
    total_area = sum(allcoral_coord$sp_area,na.rm=T)
  )
  
  for(y in years){
    y_chr <- as.character(y)
    
    # Get current year's climate values
    year_climate <- climate_dt[, .(
      x = x,
      y = y,
      dhw_max = get(paste0("dhw_max_", y_chr)),
      dhw_duration = get(paste0("dhw_duration_", y_chr))
    )]
    
    # Merge coordinates with climate data
    allcoral_pixels <- merge(allcoral_coord, year_climate, by = c("x", "y"))
    
    # Calculate exceedances
    allcoral_pixels[, `:=`(
      dhw_exceed = dhw_max > 8,
      duration_exceed = dhw_duration > 3
    )]
    
    # Calculate summary stats
    year_idx <- which(years == y)
    allcoral_summary[year_idx, `:=`(
      mean_dhw_max = weighted.mean(allcoral_pixels$dhw_max, 
                                   w = allcoral_pixels$sp_area, 
                                   na.rm = TRUE),
      mean_dhw_duration = weighted.mean(allcoral_pixels$dhw_duration, 
                                        w = allcoral_pixels$sp_area, 
                                        na.rm = TRUE),
      dhw_exceed_area = sum(allcoral_pixels[dhw_exceed == TRUE, sp_area], 
                            na.rm = TRUE),
      duration_exceed_area = sum(allcoral_pixels[duration_exceed == TRUE, sp_area], 
                                 na.rm = TRUE)
    )]
  }
  
  # Save all coral summary
  fwrite(allcoral_summary, 
         file = paste0(output_directory, "allcoral_dhwsummary_annual_", scenario, ".csv"))
}

