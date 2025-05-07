#projection analysis - threshold exceedance 

source("functions.R")

#calculate annual summer means
scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
years <- c(2020,2030,2040,2050,2060,2070,2080,2090,2100)
variables <- c("tos","ph")
models <- c('modelmean','modelmedian',"modelmeanfoc")

for(variable in variables){
  for(scenario in scenarios){
    for(model in models){
      summerannual(variable,scenario,model)
    }
  }
}






#run analysis on summer means
for(scenario in scenarios){
  #open projected climate data
  tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_lat_mask <- init(tos_summerN, fun= "y")
  tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
  
  ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_lat_mask <- init(ph_summerN, fun= "y")
  ph_summer <- ifel(ph_lat_mask >= 0, ph_summerN, ph_summerS)
  
  
  # harmonize coverage of climate data: exclude pixels where one variable is NA
  global(tos_summer, "notNA")[[1]] #43976
  global(tos_summer, "isNA")[[1]] #20824
  global(ph_summer, "notNA")[[1]] #43597 
  global(ph_summer, "isNA")[[1]] #21203
  tos_summer[is.na(ph_summer)] <- NA
  ph_summer[is.na(tos_summer)] <- NA
  #after: 43158 notNA; 21642 isNA
  
  # load species range data
  species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
  species_coord <- fread(file.path(output_directory, 'species_pixels_1982-1992.csv'))
  species_coord$x <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
  setkey(species_coord, x, y)
  
  # get unique coordinates
  coords_dt <- unique(species_coord[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))
  
  # calculate pixel areas
  area_raster <- cellSize(tos_summer, unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  # merge with species data
  species_coord <- merge(species_coord, coords_dt, by = c("x", "y"))
  species_coord[, sp_area := pixel_area * pixel_cover]
  setkey(species_coord, species_id)
  
  # extract values for each point
  tos_summer_values <- as.data.table(terra::extract(tos_summer, points, ID=F))
  ph_summer_values <- as.data.table(terra::extract(ph_summer, points, ID=F))
  
  years <- 2015:2100
  names(tos_summer_values) <- as.character(years)
  names(ph_summer_values) <- as.character(years)
  
  climate_dt <- copy(coords_dt) 
  for(i in 1:length(years)) {
    y <- years[i]
    climate_dt[, paste0("tos_", y) := tos_summer_values[[i]]]
    climate_dt[, paste0("ph_", y) := ph_summer_values[[i]]]
  }
  
  all_species_summary <-  data.table()
  
  for(sp in unique(species_envelope$id_no)){
    sp_envelope <- species_envelope[species_envelope$id_no == sp,]
    
    sp_summary <- data.table(
      species_id = sp,
      scenario = scenario,
      year = years,
      mean_temp_exceed = numeric(length(years)),  # mean of only exceeding values
      mean_ph_exceed = numeric(length(years)),       # mean of only exceeding values
      mean_temp_diff = numeric(length(years)),        # mean difference from envelope
      mean_ph_diff = numeric(length(years)),          # mean difference from envelope
      temp_exceed_area = numeric(length(years)),
      ph_exceed_area = numeric(length(years)),
      ph_exceed1_area = numeric(length(years)),       # pH more than 0.1 lower than min envelope
      ph_exceed2_area = numeric(length(years)),       # pH more than 0.2 lower than min envelope
      ph_exceed3_area = numeric(length(years)),       # pH more than 0.3 lower than min envelope
      dual_exceed_area = numeric(length(years)),
      total_area = sum(species_coord[species_id == sp, sp_area])
    )
    
    #get coords for just this species
    sp_coords <- species_coord[species_id == sp]
    
    sp_pixels_all_years <- data.table()
    
    for(y in years){
      y_chr <- as.character(y)
      
      #calculate exceedance of each pixel compared to envelope
      
      #get current year's climate values
      year_climate <- climate_dt[, .(
        x = x,
        y = y,
        tos = get(paste0("tos_", y_chr)),
        ph = get(paste0("ph_", y_chr))
      )]
      
      # get species coordinates and current year's climate values
      sp_pixels <- merge(sp_coords, year_climate, by = c("x", "y"))
      
      #calculate differences compared to envelope
      sp_pixels[, `:=`(
        year = as.numeric(y_chr),
        scenario = scenario,
        temp_diff = tos - sp_envelope$sst_max,      # all differences
        ph_diff = ph - sp_envelope$ph_min,          # all differences
        temp_exceed = pmax(0, tos - sp_envelope$sst_max),  # only exceeding
        ph_exceed = pmin(0, ph - sp_envelope$ph_min),         # only exceeding
        ph_exceed1 = pmin(0, ph - sp_envelope$ph_min + 0.1),         # exceeding 0.1
        ph_exceed2 = pmin(0, ph - sp_envelope$ph_min + 0.2),         # exceeding 0.2
        ph_exceed3 = pmin(0, ph - sp_envelope$ph_min + 0.3),         # exceeding 0.3
        dual_exceed = FALSE
      )]
      sp_pixels[temp_exceed > 0 & ph_exceed < 0, dual_exceed := TRUE]
      
      sp_pixels_all_years <- rbindlist(list(sp_pixels_all_years, 
                                            sp_pixels[, .(x, y, year, scenario, 
                                                          species_id = sp,
                                                          temp_diff,
                                                          ph_diff,
                                                          temp_exceed, 
                                                          ph_exceed, 
                                                          ph_exceed1,
                                                          ph_exceed2,
                                                          ph_exceed3,
                                                          dual_exceed,
                                                          pixel_area, 
                                                          pixel_cover,
                                                          sp_area)]))
      
      # calculate summary stats for this year
      #calculate mean exceedance per species
      #calculate total area exceeded per species
      year_idx <- which(years == y)
      sp_summary[year_idx, `:=`(
        # mean exceedance (only for exceeding pixels)
        mean_temp_exceed = weighted.mean(sp_pixels[temp_exceed > 0, temp_exceed], 
                                         w = sp_pixels[temp_exceed > 0, sp_area], 
                                         na.rm = TRUE),
        mean_ph_exceed = weighted.mean(sp_pixels[ph_exceed < 0, ph_exceed], 
                                       w = sp_pixels[ph_exceed < 0, sp_area], 
                                       na.rm = TRUE),
        # mean difference from envelope (all pixels)
        mean_temp_diff = weighted.mean(sp_pixels$temp_diff, 
                                       w = sp_pixels$sp_area, 
                                       na.rm = TRUE),
        mean_ph_diff = weighted.mean(sp_pixels$ph_diff, 
                                     w = sp_pixels$sp_area, 
                                     na.rm = TRUE),
        # areas of exceedance
        temp_exceed_area = sum(sp_pixels[temp_exceed > 0, sp_area], 
                               na.rm = TRUE),
        ph_exceed_area = sum(sp_pixels[ph_exceed < 0, sp_area], 
                             na.rm = TRUE),
        ph_exceed1_area = sum(sp_pixels[ph_exceed1 < 0, sp_area], 
                             na.rm = TRUE),
        ph_exceed2_area = sum(sp_pixels[ph_exceed2 < 0, sp_area], 
                              na.rm = TRUE),
        ph_exceed3_area = sum(sp_pixels[ph_exceed3 < 0, sp_area], 
                              na.rm = TRUE),
        dual_exceed_area = sum(sp_pixels[dual_exceed == TRUE, sp_area], 
                               na.rm = TRUE),
        total_area = sum(sp_pixels$sp_area, 
                         na.rm = TRUE)
      )]
      
    }
    
    
    # save pixel-level results for this species
    fwrite(sp_pixels_all_years, 
           file = paste0(output_directory, "threshold/species_", sp, "_exceedance_", scenario, ".csv"))
    
    # add to summary
    all_species_summary <- rbindlist(list(all_species_summary, sp_summary),fill=T)
    
  }
  
  #save summary file
  fwrite(all_species_summary, file = paste0(output_directory, "species_exceedance_summary_", scenario, ".csv"))
  
}







#run for all coral 

for(scenario in scenarios){
  allcoral_envelope <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
  allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))
  allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)
  
  
  # get unique coordinates
  coords_dt <- unique(allcoral_coord[, .(x, y)])
  setkey(coords_dt, x, y)
  points <- vect(coords_dt, geom = c("x", "y"))

  #open projected climate data
  tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_lat_mask <- init(tos_summerN, fun= "y")
  tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
  
  ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_lat_mask <- init(ph_summerN, fun= "y")
  ph_summer <- ifel(ph_lat_mask >= 0, ph_summerN, ph_summerS)
  
  
  # harmonize coverage of climate data: exclude pixels where one variable is NA
  tos_summer[is.na(ph_summer)] <- NA
  ph_summer[is.na(tos_summer)] <- NA
  
  # extract values for each point
  tos_summer_values <- as.data.table(terra::extract(tos_summer, points, ID=F))
  ph_summer_values <- as.data.table(terra::extract(ph_summer, points, ID=F))
  
  years <- 2015:2100
  names(tos_summer_values) <- as.character(years)
  names(ph_summer_values) <- as.character(years)
  
  # calculate pixel areas
  area_raster <- cellSize(tos_summer, unit='km')
  area_values <- as.data.table(terra::extract(area_raster, points))
  coords_dt[, pixel_area := area_values$area]
  
  allcoral_coord <- merge(allcoral_coord, coords_dt, by = c("x", "y"))
  allcoral_coord[, sp_area := pixel_area * pixel_cover]
  
  climate_dt <- copy(coords_dt) 
  for(i in 1:length(years)) {
    y <- years[i]
    climate_dt[, paste0("tos_", y) := tos_summer_values[[i]]]
    climate_dt[, paste0("ph_", y) := ph_summer_values[[i]]]
  }
  
  allcoral_summary <- data.table(
    species_id = "allcoral",
    scenario = scenario,
    year = years,
    mean_temp_exceed = numeric(length(years)),
    mean_ph_exceed = numeric(length(years)),
    mean_temp_diff = numeric(length(years)),
    mean_ph_diff = numeric(length(years)),
    temp_exceed_area = numeric(length(years)),
    ph_exceed_area = numeric(length(years)),
    ph_exceed1_area = numeric(length(years)),       # pH more than 0.1 lower than min envelope
    ph_exceed2_area = numeric(length(years)),       # pH more than 0.2 lower than min envelope
    ph_exceed3_area = numeric(length(years)),       # pH more than 0.3 lower than min envelope
    dual_exceed_area = numeric(length(years)),
    total_area = sum(allcoral_coord$pixel_area)
  )
  
  allcoral_pixels_all_years <- data.table()
  
  for(y in years){
    y_chr <- as.character(y)
    
    #get current year's climate values
    year_climate <- climate_dt[, .(
      x = x,
      y = y,
      tos = get(paste0("tos_", y_chr)),
      ph = get(paste0("ph_", y_chr))
    )]
    
    # merge coordinates with climate data
    allcoral_pixels <- merge(allcoral_coord, year_climate, by = c("x", "y"))
    
    # calculate exceedances
    allcoral_pixels[, `:=`(
      year = as.numeric(y_chr),
      scenario = scenario,
      temp_diff = tos - allcoral_envelope$sst_max,
      ph_diff = ph - allcoral_envelope$ph_min,
      temp_exceed = pmax(0, tos - allcoral_envelope$sst_max),
      ph_exceed = pmin(0, ph - allcoral_envelope$ph_min),
      ph_exceed1 = pmin(0, ph - allcoral_envelope$ph_min + 0.1),         # exceeding 0.1
      ph_exceed2 = pmin(0, ph - allcoral_envelope$ph_min + 0.2),         # exceeding 0.2
      ph_exceed3 = pmin(0, ph - allcoral_envelope$ph_min + 0.3),         # exceeding 0.3
      dual_exceed = FALSE
    )]
    allcoral_pixels[temp_exceed > 0 & ph_exceed < 0, dual_exceed := TRUE]
    
    # Store pixel-level data
    allcoral_pixels_all_years <- rbindlist(list(allcoral_pixels_all_years, 
                                                allcoral_pixels[, .(x, y, year, scenario, 
                                                        species_id = "allcoral",
                                                        temp_diff, ph_diff,
                                                        temp_exceed, ph_exceed, 
                                                        ph_exceed1,
                                                        ph_exceed2,
                                                        ph_exceed3,
                                                        dual_exceed,
                                                        pixel_area, 
                                                        pixel_cover,
                                                        sp_area)]))
    
    # Calculate summary statistics
    year_idx <- which(years == y)
    allcoral_summary[year_idx, `:=`(
      mean_temp_exceed = weighted.mean(allcoral_pixels[temp_exceed > 0, temp_exceed], 
                                       w = allcoral_pixels[temp_exceed > 0, sp_area], 
                                       na.rm = TRUE),
      mean_ph_exceed = weighted.mean(allcoral_pixels[ph_exceed < 0, ph_exceed], 
                                     w = allcoral_pixels[ph_exceed < 0, sp_area], 
                                     na.rm = TRUE),
      mean_temp_diff = weighted.mean(allcoral_pixels$temp_diff, 
                                     w = allcoral_pixels$sp_area, 
                                     na.rm = TRUE),
      mean_ph_diff = weighted.mean(allcoral_pixels$ph_diff, 
                                   w = allcoral_pixels$sp_area, 
                                   na.rm = TRUE),
      temp_exceed_area = sum(allcoral_pixels[temp_exceed > 0, sp_area], 
                             na.rm = TRUE),
      ph_exceed_area = sum(allcoral_pixels[ph_exceed < 0, sp_area], 
                           na.rm = TRUE),
      ph_exceed1_area = sum(allcoral_pixels[ph_exceed1 < 0, sp_area], 
                            na.rm = TRUE),
      ph_exceed2_area = sum(allcoral_pixels[ph_exceed2 < 0, sp_area], 
                            na.rm = TRUE),
      ph_exceed3_area = sum(allcoral_pixels[ph_exceed3 < 0, sp_area], 
                            na.rm = TRUE),
      dual_exceed_area = sum(allcoral_pixels[dual_exceed == TRUE, sp_area], 
                             na.rm = TRUE),
      total_area = sum(allcoral_pixels$sp_area, na.rm = TRUE)
    )]
  }
  
  # Save all coral results
  fwrite(allcoral_pixels_all_years, 
         file = paste0(output_directory, "threshold/allcoral_exceedance_", scenario, ".csv"))
  
  #save summary file
  fwrite(allcoral_summary, file = paste0(output_directory, "allcoral_exceedance_summary_", scenario, ".csv"))
}

















# realms
realm_envelopes <- fread(paste0(output_directory, 'realm_envelopes_1982-1992.csv'))
allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_with_ecoregions_1982-1992.csv'))
allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)

realms <- unique(na.omit(allcoral_coord$realm))

for(scenario in scenarios) {
  # Process each realm separately
  for(current_realm in realms) {
    # Get the envelope for the current realm
    realm_envelope <- realm_envelopes[realm == current_realm]
    
    # Filter coordinates for current realm
    realm_coord <- allcoral_coord[realm == current_realm]
    
    # Skip if no data for this realm
    if(nrow(realm_coord) == 0) next
    
    # Get unique coordinates for this realm
    coords_dt <- unique(realm_coord[, .(x, y, pixel_area, pixel_cover, sp_area)])
    setkey(coords_dt, x, y)
    points <- vect(coords_dt[, .(x, y)], geom = c("x", "y"))
    
    # Open projected climate data
    tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
    tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
    tos_lat_mask <- init(tos_summerN, fun= "y")
    tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
    
    ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
    ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
    ph_lat_mask <- init(ph_summerN, fun= "y")
    ph_summer <- ifel(ph_lat_mask >= 0, ph_summerN, ph_summerS)
    
    # Harmonize coverage of climate data: exclude pixels where one variable is NA
    tos_summer[is.na(ph_summer)] <- NA
    ph_summer[is.na(tos_summer)] <- NA
    
    # Extract values for each point
    tos_summer_values <- as.data.table(terra::extract(tos_summer, points, ID=F))
    ph_summer_values <- as.data.table(terra::extract(ph_summer, points, ID=F))
    
    years <- 2015:2100
    names(tos_summer_values) <- as.character(years)
    names(ph_summer_values) <- as.character(years)
    
    climate_dt <- copy(coords_dt) 
    for(i in 1:length(years)) {
      y <- years[i]
      climate_dt[, paste0("tos_", y) := tos_summer_values[[i]]]
      climate_dt[, paste0("ph_", y) := ph_summer_values[[i]]]
    }
    
    realm_summary <- data.table(
      species_id = "allcoral",
      realm = current_realm,
      scenario = scenario,
      year = years,
      mean_temp_exceed = numeric(length(years)),
      mean_ph_exceed = numeric(length(years)),
      mean_temp_diff = numeric(length(years)),
      mean_ph_diff = numeric(length(years)),
      temp_exceed_area = numeric(length(years)),
      ph_exceed_area = numeric(length(years)),
      ph_exceed1_area = numeric(length(years)),
      ph_exceed2_area = numeric(length(years)),
      ph_exceed3_area = numeric(length(years)),
      dual_exceed_area = numeric(length(years)),
      total_area = sum(coords_dt$sp_area, na.rm = TRUE)
    )
    
    realm_pixels_all_years <- data.table()
    
    for(y in years) {
      y_chr <- as.character(y)
      
      # Get current year's climate values
      year_climate <- climate_dt[, .(
        x = x,
        y = y,
        pixel_area = pixel_area,
        pixel_cover = pixel_cover,
        sp_area = sp_area,
        tos = get(paste0("tos_", y_chr)),
        ph = get(paste0("ph_", y_chr))
      )]
      
      # Calculate exceedances using realm-specific envelope
      year_climate[, `:=`(
        year = as.numeric(y_chr),
        scenario = scenario,
        species_id = "allcoral",
        realm = current_realm,
        temp_diff = tos - realm_envelope$sst_max,
        ph_diff = ph - realm_envelope$ph_min,
        temp_exceed = pmax(0, tos - realm_envelope$sst_max),
        ph_exceed = pmin(0, ph - realm_envelope$ph_min),
        ph_exceed1 = pmin(0, ph - realm_envelope$ph_min + 0.1),
        ph_exceed2 = pmin(0, ph - realm_envelope$ph_min + 0.2),
        ph_exceed3 = pmin(0, ph - realm_envelope$ph_min + 0.3),
        dual_exceed = FALSE
      )]
      year_climate[temp_exceed > 0 & ph_exceed < 0, dual_exceed := TRUE]
      
      # Store pixel-level data
      realm_pixels_all_years <- rbindlist(list(realm_pixels_all_years, 
                                               year_climate[, .(x, y, year, scenario, 
                                                                species_id, realm,
                                                                temp_diff, ph_diff,
                                                                temp_exceed, ph_exceed, 
                                                                ph_exceed1, ph_exceed2, ph_exceed3,
                                                                dual_exceed,
                                                                pixel_area, pixel_cover, sp_area)]))
      
      # Calculate summary statistics
      year_idx <- which(years == y)
      realm_summary[year_idx, `:=`(
        mean_temp_exceed = weighted.mean(year_climate[temp_exceed > 0, temp_exceed], 
                                         w = year_climate[temp_exceed > 0, sp_area], 
                                         na.rm = TRUE),
        mean_ph_exceed = weighted.mean(year_climate[ph_exceed < 0, ph_exceed], 
                                       w = year_climate[ph_exceed < 0, sp_area], 
                                       na.rm = TRUE),
        mean_temp_diff = weighted.mean(year_climate$temp_diff, 
                                       w = year_climate$sp_area, 
                                       na.rm = TRUE),
        mean_ph_diff = weighted.mean(year_climate$ph_diff, 
                                     w = year_climate$sp_area, 
                                     na.rm = TRUE),
        temp_exceed_area = sum(year_climate[temp_exceed > 0, sp_area], 
                               na.rm = TRUE),
        ph_exceed_area = sum(year_climate[ph_exceed < 0, sp_area], 
                             na.rm = TRUE),
        ph_exceed1_area = sum(year_climate[ph_exceed1 < 0, sp_area], 
                              na.rm = TRUE),
        ph_exceed2_area = sum(year_climate[ph_exceed2 < 0, sp_area], 
                              na.rm = TRUE),
        ph_exceed3_area = sum(year_climate[ph_exceed3 < 0, sp_area], 
                              na.rm = TRUE),
        dual_exceed_area = sum(year_climate[dual_exceed == TRUE, sp_area], 
                               na.rm = TRUE),
        total_area = sum(year_climate$sp_area, na.rm = TRUE)
      )]
    }
    
    # Save realm-specific results
    safe_realm_name <- gsub("[^a-zA-Z0-9]", "_", current_realm)  # Create safe filename
    
    fwrite(realm_pixels_all_years, 
           file = paste0(output_directory, "threshold/realm_", safe_realm_name, "_exceedance_", scenario, ".csv"))
    
    # Save summary file
    fwrite(realm_summary, 
           file = paste0(output_directory, "realm_", safe_realm_name, "_exceedance_summary_", scenario, ".csv"))
  }
  
  # Combine all realm summaries for this scenario
  combined_realm_summaries <- NULL
  for(current_realm in realms) {
    safe_realm_name <- gsub("[^a-zA-Z0-9]", "_", current_realm)
    summary_file <- paste0(output_directory, "realm_", safe_realm_name, "_exceedance_summary_", scenario, ".csv")
    
    if(file.exists(summary_file)) {
      realm_summary <- fread(summary_file)
      if(is.null(combined_realm_summaries)) {
        combined_realm_summaries <- realm_summary
      } else {
        combined_realm_summaries <- rbindlist(list(combined_realm_summaries, realm_summary))
      }
    }
  }
  
  if(!is.null(combined_realm_summaries)) {
    fwrite(combined_realm_summaries, 
           file = paste0(output_directory, "all_realms_exceedance_summary_", scenario, ".csv"))
  }
}










# ocean groups
ocean_envelopes <- fread(paste0(output_directory, 'ocean_envelopes_1982-1992.csv'))
allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_with_oceans_1982-1992.csv'))
allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)

# Get the list of oceans
oceans <- unique(na.omit(allcoral_coord$ocean))

#add indo-pacific
oceans <- c(oceans, "indo-pacific")
allcoral_coord[ocean %in% c("Indian", "Pacific"), ocean_indo_pacific := "indo-pacific"]


for(scenario in scenarios) {
  # Process each ocean separately
  for(current_ocean in oceans) {
    
    # Special handling for indo-pacific
    if(current_ocean == "indo-pacific") {
      ocean_envelope <- ocean_envelopes[ocean == current_ocean]
      ocean_coord <- allcoral_coord[ocean %in% c("Indian", "Pacific")]
    } else {
      # Original handling for other oceans
      ocean_envelope <- ocean_envelopes[ocean == current_ocean]
      ocean_coord <- allcoral_coord[ocean == current_ocean]
    }
    
    # Get unique coordinates for this ocean
    coords_dt <- unique(ocean_coord[, .(x, y, pixel_area, pixel_cover, sp_area)])
    setkey(coords_dt, x, y)
    points <- vect(coords_dt[, .(x, y)], geom = c("x", "y"))
    
    # Open projected climate data
    tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
    tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
    tos_lat_mask <- init(tos_summerN, fun= "y")
    tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
    
    ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
    ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
    ph_lat_mask <- init(ph_summerN, fun= "y")
    ph_summer <- ifel(ph_lat_mask >= 0, ph_summerN, ph_summerS)
    
    # Harmonize coverage of climate data: exclude pixels where one variable is NA
    tos_summer[is.na(ph_summer)] <- NA
    ph_summer[is.na(tos_summer)] <- NA
    
    # calculate pixel areas
    area_raster <- cellSize(tos_summer, unit='km')
    area_values <- as.data.table(terra::extract(area_raster, points))
    coords_dt[, pixel_area := area_values$area]
    coords_dt[, sp_area := pixel_area * pixel_cover]
    
    # merge with species data
    ocean_coord$pixel_area <- NULL
    ocean_coord <- merge(ocean_coord, coords_dt%>% select(x,y,pixel_area), by = c("x", "y"))
    ocean_coord[, sp_area := pixel_area * pixel_cover]

    # Extract values for each point
    tos_summer_values <- as.data.table(terra::extract(tos_summer, points, ID=F))
    ph_summer_values <- as.data.table(terra::extract(ph_summer, points, ID=F))
    
    years <- 2015:2100
    names(tos_summer_values) <- as.character(years)
    names(ph_summer_values) <- as.character(years)
    
    climate_dt <- copy(coords_dt) 
    for(i in 1:length(years)) {
      y <- years[i]
      climate_dt[, paste0("tos_", y) := tos_summer_values[[i]]]
      climate_dt[, paste0("ph_", y) := ph_summer_values[[i]]]
    }
    
    ocean_summary <- data.table(
      species_id = "allcoral",
      ocean = current_ocean,
      scenario = scenario,
      year = years,
      mean_temp_exceed = numeric(length(years)),
      mean_ph_exceed = numeric(length(years)),
      mean_temp_diff = numeric(length(years)),
      mean_ph_diff = numeric(length(years)),
      temp_exceed_area = numeric(length(years)),
      ph_exceed_area = numeric(length(years)),
      ph_exceed1_area = numeric(length(years)),
      ph_exceed2_area = numeric(length(years)),
      ph_exceed3_area = numeric(length(years)),
      dual_exceed_area = numeric(length(years)),
      total_area = sum(coords_dt$sp_area, na.rm = TRUE)
    )
    
    ocean_pixels_all_years <- data.table()
    
    for(y in years) {
      y_chr <- as.character(y)
      
      # Get current year's climate values
      year_climate <- climate_dt[, .(
        x = x,
        y = y,
        pixel_area = pixel_area,
        pixel_cover = pixel_cover,
        sp_area = sp_area,
        tos = get(paste0("tos_", y_chr)),
        ph = get(paste0("ph_", y_chr))
      )]
      
      # Calculate exceedances using ocean-specific envelope
      year_climate[, `:=`(
        year = as.numeric(y_chr),
        scenario = scenario,
        species_id = "allcoral",
        ocean = current_ocean,
        temp_diff = tos - ocean_envelope$sst_max,
        ph_diff = ph - ocean_envelope$ph_min,
        temp_exceed = pmax(0, tos - ocean_envelope$sst_max),
        ph_exceed = pmin(0, ph - ocean_envelope$ph_min),
        ph_exceed1 = pmin(0, ph - ocean_envelope$ph_min + 0.1),
        ph_exceed2 = pmin(0, ph - ocean_envelope$ph_min + 0.2),
        ph_exceed3 = pmin(0, ph - ocean_envelope$ph_min + 0.3),
        dual_exceed = FALSE
      )]
      year_climate[temp_exceed > 0 & ph_exceed < 0, dual_exceed := TRUE]
      
      # Store pixel-level data
      ocean_pixels_all_years <- rbindlist(list(ocean_pixels_all_years, 
                                               year_climate[, .(x, y, year, scenario, 
                                                                species_id, ocean,
                                                                temp_diff, ph_diff,
                                                                temp_exceed, ph_exceed, 
                                                                ph_exceed1, ph_exceed2, ph_exceed3,
                                                                dual_exceed,
                                                                pixel_area, pixel_cover, sp_area)]))
      
      # Calculate summary statistics
      year_idx <- which(years == y)
      ocean_summary[year_idx, `:=`(
        mean_temp_exceed = weighted.mean(year_climate[temp_exceed > 0, temp_exceed], 
                                         w = year_climate[temp_exceed > 0, sp_area], 
                                         na.rm = TRUE),
        mean_ph_exceed = weighted.mean(year_climate[ph_exceed < 0, ph_exceed], 
                                       w = year_climate[ph_exceed < 0, sp_area], 
                                       na.rm = TRUE),
        mean_temp_diff = weighted.mean(year_climate$temp_diff, 
                                       w = year_climate$sp_area, 
                                       na.rm = TRUE),
        mean_ph_diff = weighted.mean(year_climate$ph_diff, 
                                     w = year_climate$sp_area, 
                                     na.rm = TRUE),
        temp_exceed_area = sum(year_climate[temp_exceed > 0, sp_area], 
                               na.rm = TRUE),
        ph_exceed_area = sum(year_climate[ph_exceed < 0, sp_area], 
                             na.rm = TRUE),
        ph_exceed1_area = sum(year_climate[ph_exceed1 < 0, sp_area], 
                              na.rm = TRUE),
        ph_exceed2_area = sum(year_climate[ph_exceed2 < 0, sp_area], 
                              na.rm = TRUE),
        ph_exceed3_area = sum(year_climate[ph_exceed3 < 0, sp_area], 
                              na.rm = TRUE),
        dual_exceed_area = sum(year_climate[dual_exceed == TRUE, sp_area], 
                               na.rm = TRUE),
        total_area = sum(year_climate$sp_area, na.rm = TRUE)
      )]
    }
    
    # Save ocean-specific results
    safe_ocean_name <- gsub("[^a-zA-Z0-9]", "_", current_ocean)  # Create safe filename
    
    fwrite(ocean_pixels_all_years, 
           file = paste0(output_directory, "threshold/ocean_", safe_ocean_name, "_exceedance_", scenario, ".csv"))
    
    # Save summary file
    fwrite(ocean_summary, 
           file = paste0(output_directory, "ocean_", safe_ocean_name, "_exceedance_summary_", scenario, ".csv"))
  }
  
  # Combine all ocean summaries for this scenario
  combined_ocean_summaries <- NULL
  for(current_ocean in oceans) {
    safe_ocean_name <- gsub("[^a-zA-Z0-9]", "_", current_ocean)
    summary_file <- paste0(output_directory, "ocean_", safe_ocean_name, "_exceedance_summary_", scenario, ".csv")
    
    if(file.exists(summary_file)) {
      ocean_summary <- fread(summary_file)
      if(is.null(combined_ocean_summaries)) {
        combined_ocean_summaries <- ocean_summary
      } else {
        combined_ocean_summaries <- rbindlist(list(combined_ocean_summaries, ocean_summary))
      }
    }
  }
  
  if(!is.null(combined_ocean_summaries)) {
    fwrite(combined_ocean_summaries, 
           file = paste0(output_directory, "all_oceans_exceedance_summary_", scenario, ".csv"))
  }
}












#calculate dual exceedance as when either temp or ph exceed envelopes

scenarios <- c("ssp126","ssp245","ssp370","ssp585")

# Get list of species IDs from existing files
species_list <- unique(fread(paste0(output_directory, "species_exceedance_summary_", scenarios[1], ".csv"))$species_id)

for(scenario in scenarios) {
  # Process each species
  species_summary_file <- paste0(output_directory, "species_exceedance_summary_", scenario, ".csv")
  species_list <- unique(fread(species_summary_file)$species_id)
  
  # Initialize an empty data table to store all species summaries
  all_species_summary <- data.table()
  
  for(sp in species_list) {
    # Load existing species exceedance data
    sp_file <- paste0(output_directory, "threshold/species_", sp, "_exceedance_", scenario, ".csv")
    if(!file.exists(sp_file)) next
    
    sp_data <- fread(sp_file)
    
    # Calculate new dual exposure (EITHER temp OR pH exceeds threshold)
    sp_data[, dual_or_exceed := (temp_exceed > 0 | ph_exceed < 0)]
    
    # Calculate summary by year
    sp_summary <- sp_data[, .(
      species_id = first(species_id),
      scenario = first(scenario),
      dual_or_exceed_area = sum(sp_area[dual_or_exceed == TRUE], na.rm = TRUE)
    ), by = year]
    
    # Add to all species summary
    all_species_summary <- rbindlist(list(all_species_summary, sp_summary), fill = TRUE)
    
    # Save updated species file
    fwrite(sp_data, sp_file)
    cat("Processed species file:", sp_file, "\n")
  }
  
  # Update the summary file with calculated areas
  if(nrow(all_species_summary) > 0) {
    summary_data <- fread(species_summary_file)
    
    # Remove dual_or_exceed_area column if it exists and isn't numeric
    if("dual_or_exceed_area" %in% names(summary_data)) {
      if(!is.numeric(summary_data$dual_or_exceed_area)) {
        summary_data[, dual_or_exceed_area := NULL]
      }
    }
    
    # Merge with the summary file
    setkeyv(summary_data, c("species_id", "year", "scenario"))
    setkeyv(all_species_summary, c("species_id", "year", "scenario"))
    
    # Update or add the dual_or_exceed_area column
    summary_data[all_species_summary, dual_or_exceed_area := i.dual_or_exceed_area]
    
    # Save the updated summary file
    fwrite(summary_data, species_summary_file)
    cat("Updated species summary file:", species_summary_file, "\n")
  }
  
  # Process allcoral
  allcoral_file <- paste0(output_directory, "threshold/allcoral_exceedance_", scenario, ".csv")
  if(file.exists(allcoral_file)) {
    allcoral_data <- fread(allcoral_file)
    
    # Calculate new dual exposure (EITHER temp OR pH exceeds threshold)
    allcoral_data[, dual_or_exceed := (temp_exceed > 0 | ph_exceed < 0)]
    
    # Calculate summary by year
    allcoral_summary <- allcoral_data[, .(
      species_id = first(species_id),
      scenario = first(scenario),
      dual_or_exceed_area = sum(sp_area[dual_or_exceed == TRUE], na.rm = TRUE)
    ), by = year]
    
    # Update the summary file
    allcoral_summary_file <- paste0(output_directory, "allcoral_exceedance_summary_", scenario, ".csv")
    allcoral_summary_data <- fread(allcoral_summary_file)
    
    # Remove dual_or_exceed_area column if it exists and isn't numeric
    if("dual_or_exceed_area" %in% names(allcoral_summary_data)) {
      if(!is.numeric(allcoral_summary_data$dual_or_exceed_area)) {
        allcoral_summary_data[, dual_or_exceed_area := NULL]
      }
    }
    
    # Merge with the summary file
    setkeyv(allcoral_summary_data, c("species_id", "year"))
    setkeyv(allcoral_summary, c("species_id", "year"))
    
    # Update or add the dual_or_exceed_area column
    allcoral_summary_data[allcoral_summary, dual_or_exceed_area := i.dual_or_exceed_area]
    
    # Save the updated files
    fwrite(allcoral_data, allcoral_file)
    fwrite(allcoral_summary_data, allcoral_summary_file)
    cat("Processed allcoral for scenario:", scenario, "\n")
  }
}










