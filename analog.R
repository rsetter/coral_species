#projection analysis - climate analogs

source("functions.R")




#open species envelope data
species_coord <- read.csv(paste0(output_directory,'species_pixels_1982-1992.csv'))
species_coord$x <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
species_envelope <- read.csv(paste0(output_directory,'climate_envelopes_1982-1992.csv'))



#run for summer means
#summer analysis tells us when coral will be bleaching all summer

mask <- rast(paste0(cmip_folder,'landseamask360l.tif'))
#fix mask - all latitudes above 84.5 are NA. make them 10
lat_values <- yFromCell(mask, 1:ncell(mask))
lat_raster <- setValues(mask, lat_values)
mask[lat_raster > 83.75 & is.na(mask)] <- 10

scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
years <- c(2020,2030,2040,2050,2060,2070,2080,2090,2100)

#run analysis
for(scenario in scenarios){
  #open projected climate data
  tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OdecsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OdecsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_lat_mask <- init(tos_summerN, fun= "y")
  tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
  
  ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OdecsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OdecsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
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
  
  for(i in 1:9){ 
    run_analog_analysis(
      species_list = species_coord,
      species_envelopes = species_envelope,
      future_temp = tos_summer[[i]], 
      future_ph = ph_summer[[i]],
      time_period = years[[i]],
      base_period = "1985",
      mask = mask,
      output_dir = output_directory,
      scenario = scenario,
      model = "modelmean"
    )
    
    notifyitsdone(
      api_token = APItoken,
      channel = "coralspecies",
      event = "Analog Analysis Done",
      description = paste0("904 species analogs calculated for year ",years[[i]]),
      icon = "🔬"
    )
  }
}



### SUMMARIZE RESULTS



library(tidyverse)

for(scenario in scenarios){
  files <- list.files(output_directory,pattern = paste0("analogdiff_modelmeanfoc_",scenario,"_1985-\\d{4}\\.csv"),full.names=T)
  filesdual <- list.files(output_directory,pattern = paste0("dualanalog_modelmeanfoc_",scenario,"_1985-\\d{4}\\.csv"),full.names=T)
  data_list <- list()
  for(i in 1:length(files)) {
    year <- as.numeric(str_extract(files[i], "1985-(\\d{4})", group = 1))
    
    df <- fread(files[i])
    dfdual <- fread(filesdual[i])
    dfdual <- dfdual[,c('species_id','px','py','fx','fy','fval_temp','fval_ph','distance','directionr','directiondeg')]
    names(dfdual) <- c('species_id','px_temp','py_temp','fx_dual','fy_dual','fvaltemp_dual','fvalph_dual','distance_dual','directionr_dual','directiondeg_dual')
    df <- merge(df,dfdual,by=c("species_id","px_temp","py_temp"))
    df[, year := year]  
    data_list[[files[i]]] <- df
  }
  all_data <- rbindlist(data_list)
  
  #open species pixel cover data and join
  pixel_coverage <- fread(file.path(output_directory, 'species_pixelcover.csv'))
  names(pixel_coverage) <- c("x", "y", "pixelcover_pct","species_id")
  all_data <- merge(all_data, 
                    pixel_coverage, 
                    by.x = c("species_id", "px_temp", "py_temp"),  # adjust column names as needed
                    by.y = c("species_id", "x", "y"),
                    all.x = TRUE)
  
  all_data[is.na(pixelcover_pct), pixelcover_pct := 1] # if coverage data is missing, set pixelcover_pct to 1
  
  # Modify area calculations to use coverage fraction
  all_data[, sp_area := area * pixelcover_pct]
  
  
  # how much future condition exceeds envelopes
  # pvalf_temp - pmax_temp. pvalf_ph - pmin_ph
  all_data$pvalfExceed_temp <- all_data$pvalf_temp - all_data$pmax_temp
  all_data$pvalfExceed_ph <- all_data$pvalf_ph - all_data$pmin_ph
  
  # indicate if analog location within current distribution
  current_dist <- unique(all_data[, .(species_id, px_temp, py_temp, px_ph, py_ph)])
  all_data[, `:=`(analogInDistrib_temp = 0L,analogInDistrib_ph = 0L)]
  for (sp in unique(all_data$species_id)) {
    sp_dist <- current_dist[species_id == sp]
    all_data[species_id == sp, 
             analogInDistrib_temp := as.integer(mapply(function(x, y) 
               any(x == sp_dist$px_temp & y == sp_dist$py_temp),
               fx_temp, fy_temp))]
    all_data[species_id == sp, 
             analogInDistrib_ph := as.integer(mapply(function(x, y) 
               any(x == sp_dist$px_ph & y == sp_dist$py_ph),
               fx_ph, fy_ph))]
  }
  

  ### summarize per species per year

  summary_data <- all_data %>%
    group_by(species_id, year) %>%
    summarise(
      ### threshold exceedance 
      # pvalf_temp > pmax_temp. pvalf_ph < pmin_ph
      # calculate how many cells (and sum up the area) exceed envelope
      totalarea = sum(sp_area),
      # area where temperature/ph exceed envelope
      areaexceed_temp = sum(sp_area[pvalfExceed_temp > 0],na.rm=T),
      areaexceed_ph = sum(sp_area[pvalfExceed_ph < 0],na.rm=T), 
      # percent of area in distribution exceed envelope for temp and ph
      pctexceed_temp = areaexceed_temp / totalarea,
      pctexceed_ph = areaexceed_ph / totalarea,
      
      
      
      ### conditions at current location
      # average threshold exceedance for total distribution area of species
      pvalf_temp = mean(pvalf_temp,na.rm=T),
      pvalf_ph = mean(pvalf_ph,na.rm=T),
      aveexceed_temp = mean(pvalfExceed_temp,na.rm=T),
      aveexceed_ph = mean(pvalfExceed_ph,na.rm=T),
      
      

      ### refugias: any sites where envelope not exceeded
      # area that remains suitable
      arearefugia_temp = sum(sp_area[pvalfExceed_temp < 0],na.rm=T),
      arearefugia_ph = sum(sp_area[pvalfExceed_ph > 0],na.rm=T), 
      arearefugia_tempph = sum(sp_area[pvalfExceed_temp < 0 & pvalfExceed_ph > 0],na.rm=T),
      # percent of area in distribution that serve as refugia for temp/ph. or both
      pctrefugia_temp = arearefugia_temp / totalarea,
      pctrefugia_ph = arearefugia_ph / totalarea,
      pctrefugia_tempph = arearefugia_tempph / totalarea,
      #average adaptation requirement to other variable at refugia
      fanalogaveph_temp = mean(fanalogph_temp[pvalfExceed_temp < 0],na.rm=T), #average ph value at temp refugia
      pvalfExceedph_temp = mean(pvalfExceed_ph[pvalfExceed_temp < 0],na.rm=T), #how much average ph at temp refugia exceeds envelope 
      fanalogavetemp_ph = mean(fanalogtemp_ph[pvalfExceed_ph > 0],na.rm=T), #average temp value at ph refugia
      pvalfExceedtemp_ph = mean(pvalfExceed_temp[pvalfExceed_ph > 0],na.rm=T), #how much average temp at ph refugia exceeds envelope
      # area & percent of analogs that are still within current distribution
      refugiawithin_temp = sum(sp_area[analogInDistrib_temp == 0],na.rm=T),
      refugiawithin_ph = sum(sp_area[analogInDistrib_ph == 0],na.rm=T), 
      refugiawithin_tempph = sum(sp_area[analogInDistrib_temp == 0 & analogInDistrib_ph == 0],na.rm=T),
      pctrefugiawithin_temp = refugiawithin_temp / totalarea,
      pctrefugiawithin_ph = refugiawithin_ph / totalarea,
      pctrefugiawithin_tempph = refugiawithin_tempph / totalarea,
      
      
      
      ### analogs
      #mean distance to analog
      avedistance_temp = mean(distance_temp,na.rm=T),
      avedistance_ph = mean(distance_ph,na.rm=T),
      aveanalog_diff = mean(analog_diff,na.rm=T),
      aveanalog_dual = mean(distance_dual,na.rm=T)
      
    )
  
  #save summary csv
  fwrite(summary_data,file=paste0(output_directory,"coralsp_summary_",scenario,"_1985-2100.csv"))
  
  
}



#add on calculation - year of threshold exceedance

for(scenario in scenarios) {
  # Read the existing summary file
  summary_data <- fread(paste0(output_directory, "coralsp_summary_", scenario, "_1985-2100.csv"))
  
  # Calculate the threshold years 
  threshold_years <- summary_data %>%
    group_by(species_id) %>%
    summarise(
      firstExceed_temp = if(any(pctexceed_temp > 0, na.rm = TRUE)) 
        min(year[pctexceed_temp > 0], na.rm = TRUE) 
      else NA,
      firstExceed_ph = if(any(pctexceed_ph > 0, na.rm = TRUE)) 
        min(year[pctexceed_ph > 0], na.rm = TRUE) 
      else NA,
      completeExceed_temp = if(any(pctexceed_temp >= 0.99, na.rm = TRUE)) 
        min(year[pctexceed_temp >= 0.99], na.rm = TRUE) 
      else NA,
      completeExceed_ph = if(any(pctexceed_ph >= 0.99, na.rm = TRUE)) 
        min(year[pctexceed_ph >= 0.99], na.rm = TRUE) 
      else NA,
      completeExceed_tempph = if(any(pctexceed_temp >= 0.99 & pctexceed_ph >= 0.99, na.rm = TRUE)) 
        min(year[pctexceed_temp >= 0.99 & pctexceed_ph >= 0.99], na.rm = TRUE) 
      else NA
    )
  
  summary_data[, c("firstExceed_temp", "firstExceed_ph", 
                   "completeExceed_temp", "completeExceed_ph",
                   "completeExceed_tempph") := NULL]
  
  # Join with the correct threshold years
  summary_data <- merge(summary_data, threshold_years, by = "species_id")
  
  # Save the corrected summary data
  fwrite(summary_data, file = paste0(output_directory, "coralsp_summary_", scenario, "_1985-2100.csv"))
}










### calculate for all coral, without species breakdown

#open species envelope data
allcoral_coord <- fread(paste0(output_directory,'allcoral_pixels_1982-1992.csv'))
allcoral_coord$x <- ifelse(allcoral_coord$x < 0, allcoral_coord$x + 360, allcoral_coord$x)
allcoral_envelope <- fread(paste0(output_directory,'allcoral_envelopes_1982-1992.csv'))


#run for both summers and annual means
#summer analysis tells us when coral will be bleaching all summer
#annual means tell us when coral will be under constant state of bleaching stress

mask <- rast(paste0(cmip_folder,'landseamask360l.tif'))
#fix mask - all latitudes above 84.5 are NA. make them 10
lat_values <- yFromCell(mask, 1:ncell(mask))
lat_raster <- setValues(mask, lat_values)
mask[lat_raster > 83.75 & is.na(mask)] <- 10

scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
years <- c(2020,2030,2040,2050,2060,2070,2080,2090,2100)

#run analysis
for(scenario in scenarios){
  #open projected climate data
  tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OdecsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OdecsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_lat_mask <- init(tos_summerN, fun= "y")
  tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)
  
  ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OdecsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OdecsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
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
  
  for(i in 1:9){
    
    
    model <- "modelmean"
    time_period <- years[[i]]
    base_period <- "1985"
    future_temp <- tos_summer[[i]]
    future_ph <- ph_summer[[i]]
    
    # Create output filenames
    output_files <- list(
      dual = file.path(output_directory, 
                       sprintf("allcoral_dualanalog_%s_%s_%s-%s.csv",
                               model, scenario, base_period, time_period)),
      tos = file.path(output_directory,
                      sprintf("allcoral_analog_tos_%s_%s_%s-%s.csv",
                              model, scenario, base_period, time_period)),
      ph = file.path(output_directory,
                     sprintf("allcoral_analog_ph_%s_%s_%s-%s.csv",
                             model, scenario, base_period, time_period)),
      diff = file.path(output_directory,
                       sprintf("allcoral_analogdiff_%s_%s_%s-%s.csv",
                               model, scenario, base_period, time_period))
    )
    
    # Run dual analog analysis
    dual_results <- dualanalog(
      sp_coords = allcoral_coord,
      sp_envelope = allcoral_envelope,
      future_temp = future_temp,
      future_ph = future_ph,
      mask = mask
    )
    fwrite(dual_results, output_files$dual)
    
    
    # Run analog calc which returns a list of sst and ph results
    tosph_results <- analogcalc(
      sp_coords = allcoral_coord,
      future_temp = future_temp,
      future_ph = future_ph,
      sp_envelope = allcoral_envelope,
      mask = mask
    )
    fwrite(tosph_results$sst, output_files$tos, row.names = FALSE)
    fwrite(tosph_results$ph, output_files$ph, row.names = FALSE)
    
    # Run analog difference
    diff_results <- analogdifference(
      tempdf = tosph_results$sst,
      phdf = tosph_results$ph,
      mask = mask,
      future_temp = future_temp,
      future_ph = future_ph
    )
    fwrite(diff_results, output_files$diff, row.names = FALSE)
    
    gc()
  }
}

### SUMMARIZE RESULTS for allcoral

library(tidyverse)

for(scenario in scenarios){
  files <- list.files(output_directory,pattern = paste0("allcoral_analogdiff_modelmeanfoc_",scenario,"_1985-\\d{4}\\.csv"),full.names=T)
  filesdual <- list.files(output_directory,pattern = paste0("allcoral_dualanalog_modelmeanfoc_",scenario,"_1985-\\d{4}\\.csv"),full.names=T)
  data_list <- list()
  for(i in 1:length(files)) {
    year <- as.numeric(str_extract(files[i], "1985-(\\d{4})", group = 1))
    
    df <- fread(files[i])
    dfdual <- fread(filesdual[i])
    dfdual <- dfdual[,c('species_id','px','py','fx','fy','fval_temp','fval_ph','distance','directionr','directiondeg')]
    names(dfdual) <- c('species_id','px_temp','py_temp','fx_dual','fy_dual','fvaltemp_dual','fvalph_dual','distance_dual','directionr_dual','directiondeg_dual')
    df <- merge(df,dfdual,by=c("species_id","px_temp","py_temp"))
    df[, year := year]  
    data_list[[files[i]]] <- df
  }
  all_data <- rbindlist(data_list)
  
  #open species pixel cover data and join
  pixel_coverage <- fread(file.path(output_directory, 'allcoral_pixels_1982-1992.csv'))
  pixel_coverage <- subset(pixel_coverage,select=c("x","y","pixel_cover","species_id"))
  names(pixel_coverage) <- c("x", "y", "pixelcover_pct","species_id")
  all_data <- merge(all_data, 
                    pixel_coverage, 
                    by.x = c("species_id", "px_temp", "py_temp"),  
                    by.y = c("species_id", "x", "y"),
                    all.x = TRUE)
  all_data[is.na(pixelcover_pct), pixelcover_pct := 1] # if coverage data is missing, set pixelcover_pct to 1
  
  # Modify area calculations to use coverage fraction
  all_data[, sp_area := area * pixelcover_pct]
  
  
  # how much future condition exceeds envelopes
  # pvalf_temp - pmax_temp. pvalf_ph - pmin_ph
  all_data$pvalfExceed_temp <- all_data$pvalf_temp - all_data$pmax_temp
  all_data$pvalfExceed_ph <- all_data$pvalf_ph - all_data$pmin_ph
  
  # indicate if analog location within current distribution
  current_dist <- unique(all_data[, .(species_id, px_temp, py_temp, px_ph, py_ph)])
  all_data[, `:=`(analogInDistrib_temp = 0L,analogInDistrib_ph = 0L)]
  for (sp in unique(all_data$species_id)) {
    sp_dist <- current_dist[species_id == sp]
    all_data[species_id == sp, 
             analogInDistrib_temp := as.integer(mapply(function(x, y) 
               any(x == sp_dist$px_temp & y == sp_dist$py_temp),
               fx_temp, fy_temp))]
    all_data[species_id == sp, 
             analogInDistrib_ph := as.integer(mapply(function(x, y) 
               any(x == sp_dist$px_ph & y == sp_dist$py_ph),
               fx_ph, fy_ph))]
  }
  
  
  ### summarize per year
  
  summary_data <- all_data %>%
    group_by(year) %>%
    summarise(
      ### threshold exceedance 
      # pvalf_temp > pmax_temp. pvalf_ph < pmin_ph
      # calculate how many cells (and sum up the area) exceed envelope
      totalarea = sum(sp_area),
      # area where temperature/ph exceed envelope
      areaexceed_temp = sum(sp_area[pvalfExceed_temp > 0],na.rm=T),
      areaexceed_ph = sum(sp_area[pvalfExceed_ph < 0],na.rm=T), 
      # percent of area in distribution exceed envelope for temp and ph
      pctexceed_temp = areaexceed_temp / totalarea,
      pctexceed_ph = areaexceed_ph / totalarea,
      
      
      
      ### conditions at current location
      # average threshold exceedance for total distribution area of species
      pvalf_temp = mean(pvalf_temp,na.rm=T),
      pvalf_ph = mean(pvalf_ph,na.rm=T),
      aveexceed_temp = mean(pvalfExceed_temp,na.rm=T),
      aveexceed_ph = mean(pvalfExceed_ph,na.rm=T),
      
      
      
      ### refugias: any sites where envelope not exceeded
      # area that remains suitable
      arearefugia_temp = sum(sp_area[pvalfExceed_temp < 0],na.rm=T),
      arearefugia_ph = sum(sp_area[pvalfExceed_ph > 0],na.rm=T), 
      arearefugia_tempph = sum(sp_area[pvalfExceed_temp < 0 & pvalfExceed_ph > 0],na.rm=T),
      # percent of area in distribution that serve as refugia for temp/ph. or both
      pctrefugia_temp = arearefugia_temp / totalarea,
      pctrefugia_ph = arearefugia_ph / totalarea,
      pctrefugia_tempph = arearefugia_tempph / totalarea,
      #average adaptation requirement to other variable at refugia
      fanalogaveph_temp = mean(fanalogph_temp[pvalfExceed_temp < 0],na.rm=T), #average ph value at temp refugia
      pvalfExceedph_temp = mean(pvalfExceed_ph[pvalfExceed_temp < 0],na.rm=T), #how much average ph at temp refugia exceeds envelope 
      fanalogavetemp_ph = mean(fanalogtemp_ph[pvalfExceed_ph > 0],na.rm=T), #average temp value at ph refugia
      pvalfExceedtemp_ph = mean(pvalfExceed_temp[pvalfExceed_ph > 0],na.rm=T), #how much average temp at ph refugia exceeds envelope
      # area & percent of analogs that are still within current distribution
      refugiawithin_temp = sum(sp_area[analogInDistrib_temp == 0],na.rm=T),
      refugiawithin_ph = sum(sp_area[analogInDistrib_ph == 0],na.rm=T), 
      refugiawithin_tempph = sum(sp_area[analogInDistrib_temp == 0 & analogInDistrib_ph == 0],na.rm=T),
      pctrefugiawithin_temp = refugiawithin_temp / totalarea,
      pctrefugiawithin_ph = refugiawithin_ph / totalarea,
      pctrefugiawithin_tempph = refugiawithin_tempph / totalarea,
      
      
      
      ### analogs
      #mean distance to analog
      avedistance_temp = mean(distance_temp,na.rm=T),
      avedistance_ph = mean(distance_ph,na.rm=T),
      aveanalog_diff = mean(analog_diff,na.rm=T),
      aveanalog_dual = mean(distance_dual,na.rm=T)
      
    )
  
  #save summary csv
  fwrite(summary_data,file=paste0(output_directory,"allcoral_summary_",scenario,"_1985-2100.csv"))
  
  
}






