library(parallel)
library(doParallel)
library(terra)
library(data.table)


# for each species, what's the min, mean, median, max for its full distribution? - for temp and ph (max = maximum monthly mean)
# maximum monthly mean (MMM) was calculated as the maximum monthly SST between 1982 and 1992
# Then when the mean summer temp exceeds the max monthly mean, that means that coral is bleaching the entire summer, or 12 weeks of DHW


# open historic cmip6 monthly multi-model mean 
sst <- rast(paste0(cmip_folder,'tos/tos_Omon_modelmean_hist_185001-201412.tif'))
ph <- rast(paste0(cmip_folder,'ph/ph_Omon_modelmean_hist_185001-201412.tif'))

#rotate 0,360 to -180,180 to match polygons
sst180 <- rotate(sst)
sst180r <- crop(sst180, ext(-180.5, -180, -90, 90))
sst180l <- crop(sst180, ext(-180, 179.5, -90, 90))
sst180r <- shift(sst180r, dx = 360)
sst180 <- merge(sst180l, sst180r)

ph180 <- rotate(ph)
ph180r <- crop(ph180, ext(-180.5, -180, -90, 90))
ph180l <- crop(ph180, ext(-180, 179.5, -90, 90))
ph180r <- shift(ph180r, dx = 360)
ph180 <- merge(ph180l, ph180r)

# calculate climate stats - min, mean, median, max for each pixel across the time period for sst and ph
# 1982-1992 period
start_8292 <- (1982-1850)*12 + 1
end_8292 <- (1992-1850)*12 + 12
sst_8292 <- sst180[[start_8292:end_8292]]
ph_8292 <- ph180[[start_8292:end_8292]]

sst_8292_min <- min(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1982-1992_min.tif'))
sst_8292_max <- max(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1982-1992_max.tif'))
sst_8292_mean <- mean(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1982-1992_mean.tif'))
sst_8292_median <- app(sst_8292, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1982-1992_median.tif'))

ph_8292_min <- min(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1982-1992_min.tif'))
ph_8292_max <- max(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1982-1992_max.tif'))
ph_8292_mean <- mean(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1982-1992_mean.tif'))
ph_8292_median <- app(ph_8292, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1982-1992_median.tif'))


# 1850-1990 period
start_hist <- 1
end_hist <- (1990-1850)*12 + 12
sst_hist <- sst180[[start_hist:end_hist]]
ph_hist <- ph180[[start_hist:end_hist]]

sst_hist_min <- min(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1850-1990_min.tif'))
sst_hist_max <- max(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1850-1990_max.tif'))
sst_hist_mean <- mean(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1850-1990_mean.tif'))
sst_hist_median <- app(sst_hist, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_1850-1990_median.tif'))

ph_hist_min <- min(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1850-1990_min.tif'))
ph_hist_max <- max(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1850-1990_max.tif'))
ph_hist_mean <- mean(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1850-1990_mean.tif'))
ph_hist_median <- app(ph_hist, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_1850-1990_median.tif'))



#for each species, find min, mean, median, max within its polygon for sst and ph

periods <- list(
  '8292' = list(
    sst_files = c(
      file.path(iucn_sp_folder, 'sst_1982-1992_min.tif'),
      file.path(iucn_sp_folder, 'sst_1982-1992_max.tif'),
      file.path(iucn_sp_folder, 'sst_1982-1992_mean.tif'),
      file.path(iucn_sp_folder, 'sst_1982-1992_median.tif')
    ),
    ph_files = c(
      file.path(iucn_sp_folder, 'ph_1982-1992_min.tif'),
      file.path(iucn_sp_folder, 'ph_1982-1992_max.tif'),
      file.path(iucn_sp_folder, 'ph_1982-1992_mean.tif'),
      file.path(iucn_sp_folder, 'ph_1982-1992_median.tif')
    ),
    outfile = file.path(output_directory, 'climate_envelopes_1982-1992.csv'),
    coords_file = file.path(output_directory, 'species_pixels_1982-1992.csv')
  ),
  'hist' = list(
    sst_files = c(
      file.path(iucn_sp_folder, 'sst_1850-1990_min.tif'),
      file.path(iucn_sp_folder, 'sst_1850-1990_max.tif'),
      file.path(iucn_sp_folder, 'sst_1850-1990_mean.tif'),
      file.path(iucn_sp_folder, 'sst_1850-1990_median.tif')
    ),
    ph_files = c(
      file.path(iucn_sp_folder, 'ph_1850-1990_min.tif'),
      file.path(iucn_sp_folder, 'ph_1850-1990_max.tif'),
      file.path(iucn_sp_folder, 'ph_1850-1990_mean.tif'),
      file.path(iucn_sp_folder, 'ph_1850-1990_median.tif')
    ),
    outfile = file.path(output_directory, 'climate_envelopes_1850-1990.csv'),
    coords_file = file.path(output_directory, 'species_pixels_1850-1990.csv')
  )
)

# Set up parallel processing
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(terra)
  library(data.table)
  TRUE
})

# Set chunk size
chunk_size <- 50
n_species <- 903

# Process each time period
for(period_name in names(periods)) {
  cat("\nProcessing", period_name, "period\n")
  period <- periods[[period_name]]
  
  # Process chunks
  for(i in seq(1, n_species, chunk_size)) {
    end_i <- min(i + chunk_size - 1, n_species)
    cat("Processing species", i, "to", end_i, "\n")
    
    # Parallel processing within the chunk
    chunk_results <- foreach(j = i:end_i, 
                             .packages = c("terra", "data.table")) %dopar% {
                               
                               #open species polygon
                               species_query <- sprintf("SELECT * FROM coral_all WHERE fid = %d", j)
                               species_poly <- vect(paste0(iucn_sp_folder, "coral_all.gpkg"), query=species_query)
                               species_id <- species_poly$id_no
                               
                               # open climate rasters
                               sst_stats <- rast(period$sst_files)
                               ph_stats <- rast(period$ph_files)
                               
                               # Extract values
                               sst_values <- terra::extract(sst_stats, species_poly,xy=T)
                               ph_values <- terra::extract(ph_stats, species_poly,xy=T)
                               
                               # Get column names for the extracted values
                               sst_cols <- names(sst_values)[!(names(sst_values) %in% c("ID", "x", "y"))]
                               ph_cols <- names(ph_values)[!(names(ph_values) %in% c("ID", "x", "y"))]
                               
                               # compile pixel-level data for each species
                               pixel_data <- data.frame(
                                 species_id = species_id,
                                 x = sst_values$x,
                                 y = sst_values$y,
                                 pixel_id = 1:nrow(sst_values),
                                 sst_min = sst_values[[sst_cols[1]]],
                                 sst_max = sst_values[[sst_cols[2]]],
                                 sst_mean = sst_values[[sst_cols[3]]],
                                 sst_median = sst_values[[sst_cols[4]]],
                                 ph_min = ph_values[[ph_cols[1]]],
                                 ph_max = ph_values[[ph_cols[2]]],
                                 ph_mean = ph_values[[ph_cols[3]]],
                                 ph_median = ph_values[[ph_cols[4]]]
                               )
                               
                               # Calculate species-level summary statistics
                               species_data <- as.data.frame(species_poly)
                               climate_stats <- data.frame(
                                 n_pixels = nrow(sst_values),
                                 sst_min = min(sst_values[[sst_cols[1]]], na.rm = TRUE),
                                 sst_max = max(sst_values[[sst_cols[2]]], na.rm = TRUE),
                                 sst_mean = mean(sst_values[[sst_cols[3]]], na.rm = TRUE),
                                 sst_median = median(sst_values[[sst_cols[4]]], na.rm = TRUE),
                                 ph_min = min(ph_values[[ph_cols[1]]], na.rm = TRUE),
                                 ph_max = max(ph_values[[ph_cols[2]]], na.rm = TRUE),
                                 ph_mean = mean(ph_values[[ph_cols[3]]], na.rm = TRUE),
                                 ph_median = median(ph_values[[ph_cols[4]]], na.rm = TRUE)
                               )
                               
                               species_stats <- cbind(species_data, climate_stats)
                               
                               list(
                                 species_stats = species_stats,
                                 pixel_data = pixel_data
                               )
                             }
    # Separate species stats and pixel data
    species_stats_df <- do.call(rbind, lapply(chunk_results, function(x) x$species_stats))
    pixel_data_df <- do.call(rbind, lapply(chunk_results, function(x) x$pixel_data))
    
    # Convert to data.table for fwrite
    species_stats_dt <- as.data.table(species_stats_df)
    pixel_data_dt <- as.data.table(pixel_data_df)
    
    # Write species-level climate envelopes
    fwrite(species_stats_dt, 
           file = period$outfile, 
           append = i != 1,  # append if not first chunk
           col.names = i == 1)  # write headers only for first chunk
    
    # Write pixel-level data
    fwrite(pixel_data_dt, 
           file = period$coords_file,
           append = i != 1,
           col.names = i == 1)
    
    rm(chunk_results, species_stats_df, pixel_data_df,species_stats_dt,pixel_data_dt)
    gc()
  }
}

stopCluster(cl)


species_coord <- read.csv(periods[['8292']]$coords_file)
species_envelope <- read.csv(periods[['8292']]$outfile)







#calculate fractional coverage for each pixel
cores <- detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(terra)
  library(data.table)
  TRUE
})

# Set chunk size
chunk_size <- 50
n_species <- 903

for(i in seq(1, n_species, chunk_size)) {
  end_i <- min(i + chunk_size - 1, n_species)
  cat("Processing species", i, "to", end_i, "\n")
  
  chunk_results <- foreach(j = i:end_i, 
                           .packages = c("terra")) %dopar% {
                             
                             # Open species polygon
                             species_query <- sprintf("SELECT * FROM coral_all WHERE fid = %d", j)
                             species_poly <- vect(paste0(iucn_sp_folder, "coral_all.gpkg"), query=species_query)
                             species_id <- species_poly$id_no
                             
                             # Create template raster covering species extent
                             template_rast <- rast(periods[['8292']]$sst_files[1])  # Using first file for resolution
                             species_rast <- rasterize(species_poly, template_rast, cover=TRUE)
                             
                             # Extract cells with coverage values
                             cells <- as.data.frame(species_rast, xy=TRUE)
                             
                             # Remove NA values and rename columns
                             pixel_coverage <- na.omit(cells)
                             names(pixel_coverage) <- c("x", "y", "coverage_fraction")
                             
                             # Add species ID
                             pixel_coverage$species_id <- species_id
                             
                             pixel_coverage
                           }
  #combine into data.table, save
  pixelcover <- rbindlist(chunk_results)
  fwrite(pixelcover,
         file= file.path(output_directory, 'species_pixelcover.csv'),
         append = i != 1,
         col.names = i == 1)
  
  #clean
  rm(chunk_results,pixelcover)
  gc()
}







library(RCurl) #notification system - phone notification when code finished running

# set up notification system
headers = c(
  "Content-Type" = "application/json",
  "Authorization" = APItoken #insert personal API token (generated when sign up)
)
params = "{
  \"project\": \"rcode\",
  \"channel\": \"rcode\",
  \"event\": \"analog-calc is complete\",
  \"description\": \"code finished\",
  \"icon\": \"🔥\",
  \"notify\": true
}"
#code below sends the notification
res <- postForm("https://api.logsnag.com/v1/log", .opts=list(postfields = params, httpheader = headers, followlocation = TRUE), style = "httppost") 








