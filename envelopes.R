library(parallel)
library(doParallel)
library(terra)
library(data.table)


source("functions.R")

# for each species, what's the min, mean, median, max for its full distribution? - for temp and ph (max = maximum monthly mean)
# maximum monthly mean (MMM) was calculated as the maximum monthly SST between 1982 and 1992
# Then when the mean summer temp exceeds the max monthly mean, that means that coral is bleaching the entire summer, or 12 weeks of DHW


# open historic cmip6 monthly multi-model mean 
sst <- rast(paste0(cmip_folder,'tos/tos_Omon_modelmeanfoc_hist_185001-201412.tif'))
ph <- rast(paste0(cmip_folder,'ph/ph_Omon_modelmeanfoc_hist_185001-201412.tif'))

# harmonize coverage of climate data: exclude pixels where one variable is NA
global(sst, "notNA")[[1]] #47506 
global(sst, "isNA")[[1]] #17294 
global(ph, "notNA")[[1]] #46978 
global(ph, "isNA")[[1]] #17822 
sst[is.na(ph)] <- NA
ph[is.na(sst)] <- NA

#rotate 0,360 to -180,180 to match polygons
sst180 <- rotate(sst)
sst180r <- crop(sst180, ext(-180.5, -180, -90, 90))
sst180l <- crop(sst180, ext(-180, 179.5, -90, 90))
sst180r <- terra::shift(sst180r, dx = 360)
sst180 <- merge(sst180l, sst180r)

ph180 <- rotate(ph)
ph180r <- crop(ph180, ext(-180.5, -180, -90, 90))
ph180l <- crop(ph180, ext(-180, 179.5, -90, 90))
ph180r <- terra::shift(ph180r, dx = 360)
ph180 <- merge(ph180l, ph180r)

# calculate climate stats - min, mean, median, max for each pixel across the time period for sst and ph
# 1982-1992 period
start_8292 <- (1982-1850)*12 + 1
end_8292 <- (1992-1850)*12 + 12
sst_8292 <- sst180[[start_8292:end_8292]]
ph_8292 <- ph180[[start_8292:end_8292]]

sst_8292_min <- min(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1982-1992_min.tif'))
sst_8292_max <- max(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1982-1992_max.tif'))
sst_8292_mean <- mean(sst_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1982-1992_mean.tif'))
sst_8292_median <- app(sst_8292, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1982-1992_median.tif'))

ph_8292_min <- min(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1982-1992_min.tif'))
ph_8292_max <- max(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1982-1992_max.tif'))
ph_8292_mean <- mean(ph_8292, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1982-1992_mean.tif'))
ph_8292_median <- app(ph_8292, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1982-1992_median.tif'))


# 1850-1990 period
start_hist <- 1
end_hist <- (1990-1850)*12 + 12
sst_hist <- sst180[[start_hist:end_hist]]
ph_hist <- ph180[[start_hist:end_hist]]

sst_hist_min <- min(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1850-1990_min.tif'))
sst_hist_max <- max(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1850-1990_max.tif'))
sst_hist_mean <- mean(sst_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1850-1990_mean.tif'))
sst_hist_median <- app(sst_hist, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'sst_foc_1850-1990_median.tif'))

ph_hist_min <- min(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1850-1990_min.tif'))
ph_hist_max <- max(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1850-1990_max.tif'))
ph_hist_mean <- mean(ph_hist, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1850-1990_mean.tif'))
ph_hist_median <- app(ph_hist, median, na.rm=TRUE,filename=paste0(iucn_sp_folder,'ph_foc_1850-1990_median.tif'))



#for each species, find min, mean, median, max within its polygon for sst and ph

periods <- list(
  '8292' = list(
    sst_files = c(
      file.path(iucn_sp_folder, 'sst_foc_1982-1992_min.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1982-1992_max.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1982-1992_mean.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1982-1992_median.tif')
    ),
    ph_files = c(
      file.path(iucn_sp_folder, 'ph_foc_1982-1992_min.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1982-1992_max.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1982-1992_mean.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1982-1992_median.tif')
    ),
    outfile = file.path(output_directory, 'climate_envelopes_1982-1992.csv'),
    coords_file = file.path(output_directory, 'species_pixels_1982-1992.csv')
  ),
  'hist' = list(
    sst_files = c(
      file.path(iucn_sp_folder, 'sst_foc_1850-1990_min.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1850-1990_max.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1850-1990_mean.tif'),
      file.path(iucn_sp_folder, 'sst_foc_1850-1990_median.tif')
    ),
    ph_files = c(
      file.path(iucn_sp_folder, 'ph_foc_1850-1990_min.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1850-1990_max.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1850-1990_mean.tif'),
      file.path(iucn_sp_folder, 'ph_foc_1850-1990_median.tif')
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
                               sst_values <- terra::extract(sst_stats, species_poly,xy=T,exact=T)
                               ph_values <- terra::extract(ph_stats, species_poly,xy=T,exact=T)
                               
                               # Get column names for the extracted values
                               sst_cols <- names(sst_values)[!(names(sst_values) %in% c("ID", "x", "y","fraction"))]
                               ph_cols <- names(ph_values)[!(names(ph_values) %in% c("ID", "x", "y","fraction"))]
                               
                               # compile pixel-level data for each species
                               pixel_data <- data.frame(
                                 species_id = species_id,
                                 x = sst_values$x,
                                 y = sst_values$y,
                                 pixel_cover = sst_values$fraction,
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
                               # use 95th percentile to exclude any outliers (2.5 for min and 97.5 for max)
                               species_data <- as.data.frame(species_poly)
                               climate_stats <- data.frame(
                                 n_pixels = nrow(sst_values),
                                 sst_min_absolute = min(sst_values[[sst_cols[1]]], na.rm = TRUE),
                                 sst_max_absolute = max(sst_values[[sst_cols[2]]], na.rm = TRUE),
                                 sst_min = quantile(sst_values[[sst_cols[1]]], probs = 0.025, na.rm = TRUE),
                                 sst_max = quantile(sst_values[[sst_cols[2]]], probs = 0.975, na.rm = TRUE),
                                 sst_mean = mean(sst_values[[sst_cols[3]]], na.rm = TRUE),
                                 sst_median = median(sst_values[[sst_cols[4]]], na.rm = TRUE),
                                 ph_min_absolute = min(ph_values[[ph_cols[1]]], na.rm = TRUE),
                                 ph_max_absolute = max(ph_values[[ph_cols[2]]], na.rm = TRUE),
                                 ph_min = quantile(ph_values[[ph_cols[1]]], probs = 0.025, na.rm = TRUE),
                                 ph_max = quantile(ph_values[[ph_cols[2]]], probs = 0.975, na.rm = TRUE),
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








#calculate overall coral envelope & distribution
species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
species_coord <- fread(paste0(output_directory, 'species_pixels_1982-1992.csv'))

allcoral_envelope <- data.frame(
  id_no = "allcoral",
  sst_min_absolute = min(species_envelope$sst_min, na.rm = TRUE),
  sst_max_absolute = max(species_envelope$sst_max, na.rm = TRUE),
  sst_min = quantile(species_envelope$sst_min, probs = 0.025, na.rm = TRUE),
  sst_max = quantile(species_envelope$sst_max, probs = 0.975, na.rm = TRUE),
  sst_mean = mean(species_envelope$sst_mean, na.rm = TRUE),
  sst_median = mean(species_envelope$sst_median, na.rm = TRUE),
  ph_min_absolute = min(species_envelope$ph_min, na.rm = TRUE),
  ph_max_absolute = max(species_envelope$ph_max, na.rm = TRUE),
  ph_min = quantile(species_envelope$ph_min, probs = 0.025, na.rm = TRUE),
  ph_max = quantile(species_envelope$ph_max, probs = 0.975, na.rm = TRUE),
  ph_mean = mean(species_envelope$ph_mean, na.rm = TRUE),
  ph_median = mean(species_envelope$ph_median, na.rm = TRUE)
)

allcoral_coord <- species_coord %>%
  group_by(x, y) %>%
  summarize(
    species_id = "allcoral",
    pixel_cover = max(pixel_cover, na.rm = TRUE),
    pixel_id = first(pixel_id),  
    sst_min = first(sst_min),   
    sst_max = first(sst_max),    
    sst_mean = first(sst_mean),
    sst_median = first(sst_median),
    ph_min = first(ph_min),
    ph_max = first(ph_max),
    ph_mean = first(ph_mean),
    ph_median = first(ph_median)
  ) %>%
  ungroup()

fwrite(allcoral_envelope,paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
fwrite(allcoral_coord, paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))






## add summer means

species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
species_coord <- fread(paste0(output_directory, 'species_pixels_1982-1992.csv'))

allcoral_envelope <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))

sst <- rast(paste0(cmip_folder,'tos/tos_Omon_modelmeanfoc_hist_185001-201412.tif'))
ph <- rast(paste0(cmip_folder,'ph/ph_Omon_modelmeanfoc_hist_185001-201412.tif'))
start_idx <- (1982-1850)*12 + 1
end_idx <- (1992-1850)*12 + 12
months_8292 <- start_idx:end_idx
sst_8292 <- sst[[months_8292]]
ph_8292 <- ph[[months_8292]]

n_years <- 11
sst_summerN_layers <- c()
sst_summerS_layers <- c()
ph_summerN_layers <- c()
ph_summerS_layers <- c()

# Northern hemisphere summer (Jul-Sep), Southern hemisphere summer (Jan-Mar)
for (year in 0:(n_years-1)) {
  base_idx <- year * 12
  
  # Northern hemisphere summer indices (Jul-Sep)
  jul_idx <- base_idx + 7
  aug_idx <- base_idx + 8
  sep_idx <- base_idx + 9
  
  # Southern hemisphere summer indices (Jan-Mar)
  jan_idx <- base_idx + 1
  feb_idx <- base_idx + 2
  mar_idx <- base_idx + 3
  
  # Collect all layers
  sst_summerN_layers <- c(sst_summerN_layers, c(jul_idx, aug_idx, sep_idx))
  sst_summerS_layers <- c(sst_summerS_layers, c(jan_idx, feb_idx, mar_idx))
  ph_summerN_layers <- c(ph_summerN_layers, c(jul_idx, aug_idx, sep_idx))
  ph_summerS_layers <- c(ph_summerS_layers, c(jan_idx, feb_idx, mar_idx))
}

sst_summerN <- sst_8292[[sst_summerN_layers]]
sst_summerS <- sst_8292[[sst_summerS_layers]]
ph_summerN <- ph_8292[[ph_summerN_layers]]
ph_summerS <- ph_8292[[ph_summerS_layers]]

sst_summerN_mean <- mean(sst_summerN, na.rm=TRUE)
sst_summerS_mean <- mean(sst_summerS, na.rm=TRUE)
ph_summerN_mean <- mean(ph_summerN, na.rm=TRUE)
ph_summerS_mean <- mean(ph_summerS, na.rm=TRUE)

sst_lat_mask <- init(sst_summerN_mean, fun= "y")
ph_lat_mask <- init(ph_summerN_mean, fun= "y")

sst_summer_mean <- ifel(sst_lat_mask >= 0, sst_summerN_mean, sst_summerS_mean)
ph_summer_mean <- ifel(ph_lat_mask >= 0, ph_summerN_mean, ph_summerS_mean)

sst180 <- rotate(sst_summer_mean)
sst180r <- crop(sst180, ext(-180.5, -180, -90, 90))
sst180l <- crop(sst180, ext(-180, 179.5, -90, 90))
sst180r <- terra::shift(sst180r, dx = 360)
sst_summer_mean <- merge(sst180l, sst180r)

ph180 <- rotate(ph_summer_mean)
ph180r <- crop(ph180, ext(-180.5, -180, -90, 90))
ph180l <- crop(ph180, ext(-180, 179.5, -90, 90))
ph180r <- terra::shift(ph180r, dx = 360)
ph_summer_mean <- merge(ph180l, ph180r)

coords <- unique(species_coord[, .(x, y)])
points <- vect(coords, geom=c("x", "y"), crs=crs(sst_summer_mean))

# Extract summer values at each point
summer_values <- terra::extract(c(sst_summer_mean, ph_summer_mean), points)

# Create a lookup table with coordinates and summer values
lookup_table <- data.table(
  x = coords$x,
  y = coords$y,
  sst_summer_mean = summer_values[[2]],
  ph_summer_mean = summer_values[[3]]
)

# Merge summer values into the species coordinates data
species_coord <- merge(species_coord, lookup_table, by=c("x", "y"), all.x=TRUE)
allcoral_coord <- merge(allcoral_coord, lookup_table, by=c("x", "y"), all.x=TRUE)

# Calculate species-level summer statistics
species_summer <- species_coord[, .(
  sst_summer_mean = mean(sst_summer_mean, na.rm=TRUE),
  ph_summer_mean = mean(ph_summer_mean, na.rm=TRUE)
), by=species_id]
allcoral_summer <- allcoral_coord[, .(
  sst_summer_mean = mean(sst_summer_mean, na.rm=TRUE),
  ph_summer_mean = mean(ph_summer_mean, na.rm=TRUE)
), by=species_id]

# Merge summer statistics into the species envelope data
species_envelope <- merge(species_envelope, species_summer,by.x="id_no", by.y="species_id", all.x=TRUE)
allcoral_envelope <- merge(allcoral_envelope, allcoral_summer,by.x="id_no", by.y="species_id", all.x=TRUE)


fwrite(species_coord,paste0(output_directory, 'species_pixels_1982-1992.csv'))
fwrite(species_envelope,paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
fwrite(allcoral_envelope,paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
fwrite(allcoral_coord, paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))
