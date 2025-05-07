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

















###  ecoregion scale 

meow <- vect(paste0(ecosystem_folder, "Marine Ecoregions of the World"))
coral_points <- vect(species_coord, geom=c("x", "y"), crs=crs(meow))
ecosystems_df <- terra::extract(meow, coral_points) #15% of pixels not assigned a realm


# Add ecosystem information to the species coordinates
species_coord$realm <- ecosystems_df$REALM
species_coord$province <- ecosystems_df$PROVINCE
species_coord$ecoregion <- ecosystems_df$ECOREGION

ecoregion_lookup <- species_coord[, .(x, y, realm, province, ecoregion)] %>% unique()
allcoral_coord <- merge(allcoral_coord, ecoregion_lookup, by=c("x", "y"), all.x=TRUE)


# Function to calculate envelopes for a group of points
calculate_envelope <- function(group_data, group_name, level_name) {
  envelope <- data.frame(
    id_no = paste0(group_name, "_", level_name),
    realm = level_name,
    sst_min_absolute = min(group_data$sst_min, na.rm = TRUE),
    sst_max_absolute = max(group_data$sst_max, na.rm = TRUE),
    sst_min = quantile(group_data$sst_min, probs = 0.025, na.rm = TRUE),
    sst_max = quantile(group_data$sst_max, probs = 0.975, na.rm = TRUE),
    sst_mean = mean(group_data$sst_mean, na.rm = TRUE),
    sst_median = mean(group_data$sst_median, na.rm = TRUE),
    sst_summer_mean = mean(group_data$sst_summer_mean, na.rm = TRUE),
    ph_min_absolute = min(group_data$ph_min, na.rm = TRUE),
    ph_max_absolute = max(group_data$ph_max, na.rm = TRUE),
    ph_min = quantile(group_data$ph_min, probs = 0.025, na.rm = TRUE),
    ph_max = quantile(group_data$ph_max, probs = 0.975, na.rm = TRUE),
    ph_mean = mean(group_data$ph_mean, na.rm = TRUE),
    ph_median = mean(group_data$ph_median, na.rm = TRUE),
    ph_summer_mean = mean(group_data$ph_summer_mean, na.rm = TRUE)
  )
  return(envelope)
}

# Calculate envelopes by realm for all corals combined
realms <- unique(allcoral_coord$realm)
realm_envelopes <- data.table()

for (realm_name in realms) {
  if (is.na(realm_name)) next
  
  realm_data <- allcoral_coord[realm == realm_name]
  if (nrow(realm_data) > 0) {
    realm_envelope <- calculate_envelope(realm_data, "allcoral", realm_name)
    realm_envelopes <- rbind(realm_envelopes, realm_envelope)
  }
}


# Write results to CSV
fwrite(realm_envelopes, paste0(output_directory, 'realm_envelopes_1982-1992.csv'))

# Update the main dataframes with realm information
fwrite(species_coord, paste0(output_directory, 'species_pixels_with_ecoregions_1982-1992.csv'))
fwrite(allcoral_coord, paste0(output_directory, 'allcoral_pixels_with_ecoregions_1982-1992.csv'))










#oceans
oceans <- vect(paste0(ecosystem_folder, "GlobalOceans"))


allcoral_coord180 <- allcoral_coord

# Convert longitude from 0-360° to -180 to 180°
allcoral_coord180$x_180 <- ifelse(allcoral_coord$x > 180, 
                                   allcoral_coord$x - 360, 
                                   allcoral_coord$x)

coral_pt <- vect(allcoral_coord180, geom=c("x_180", "y"), crs="EPSG:4326")

coral_oceans <- extract(oceans, coral_pt)

# Join the ocean basin information back to original data frame
allcoral_coord$ocean_basin <- coral_oceans$name
sum(is.na(allcoral_coord$ocean_basin))
na_rows <- which(is.na(allcoral_coord$ocean_basin))

# For each NA point, find the closest point with a non-NA ocean_basin
for (i in na_rows) {
  x_i <- allcoral_coord$x[i]
  y_i <- allcoral_coord$y[i]
  
  # Define what "adjacent" means - points within 1 degree
  adjacent_threshold <- 1
  
  # Find non-NA points
  non_na_points <- which(!is.na(allcoral_coord$ocean_basin))
  
  # Calculate distances to all non-NA points
  distances <- sqrt((allcoral_coord$x[non_na_points] - x_i)^2 + 
                      (allcoral_coord$y[non_na_points] - y_i)^2)
  
  # Find the closest point within threshold
  adjacent_indices <- non_na_points[distances <= adjacent_threshold]
  
  if (length(adjacent_indices) > 0) {
    # Get the most common ocean_basin among adjacent points
    adjacent_groups <- allcoral_coord$ocean_basin[adjacent_indices]
    realm_group_counts <- table(adjacent_groups)
    most_common_group <- names(realm_group_counts)[which.max(realm_group_counts)]
    
    # Assign the most common adjacent ocean_basin
    allcoral_coord$ocean_basin[i] <- most_common_group
  }
}


unique(allcoral_coord$ocean_basin)

allcoral_coord$ocean <- allcoral_coord$ocean_basin
allcoral_coord$ocean <- gsub("North Pacific.*|South Pacific.*|.*Pacific.*", "Pacific", allcoral_coord$ocean)
allcoral_coord$ocean <- gsub(".*South China and Easter Archipelagic Seas.*", "Pacific", allcoral_coord$ocean)
allcoral_coord$ocean <- gsub(".*Mediterranean Region.*" , "Atlantic", allcoral_coord$ocean)
allcoral_coord$ocean <- gsub("North Atlantic.*|South Atlantic.*|.*Atlantic.*", "Atlantic", allcoral_coord$ocean)
allcoral_coord$ocean <- gsub(".*Indian.*", "Indian", allcoral_coord$ocean)

# fix ocean-pixel assignmetns due to mismatch in resolution
species_per_pixel <- species_coord[, .(num_species = uniqueN(species_id)), by = .(x, y)]
ocean_lookup <- allcoral_coord[, .(lon, y, current_ocean = ocean)]
names(ocean_lookup) <- c("x","y","current_ocean")
species_locations <- unique(species_coord[, .(x, y, species_id)])
species_with_ocean <- merge(species_locations, ocean_lookup, by = c("x", "y"), all.x = TRUE)
species_ocean_dist <- species_with_ocean[!is.na(current_ocean), .(pixel_count = .N), by = .(species_id, current_ocean)]
species_totals <- species_ocean_dist[, .(total_pixels = sum(pixel_count)), by = species_id]
species_ocean_pct <- merge(species_ocean_dist, species_totals, by = "species_id")
species_ocean_pct[, pct_in_ocean := pixel_count / total_pixels]
species_ocean_wide <- dcast(species_ocean_pct, species_id ~ current_ocean,value.var = "pct_in_ocean", fill = 0)
setnames(species_ocean_wide,  c("Pacific", "Atlantic", "Indian"), c("pct_range_pac", "pct_range_atl", "pct_range_ind"), skip_absent = TRUE)
species_per_pixel <- species_locations[, .(num_species = .N), by = .(x, y)]
pixel_species <- species_locations[, .(species_list = list(species_id)), by = .(x, y)]

calculate_ocean_pcts <- function(species_list, species_ocean_wide) {
  subset_data <- species_ocean_wide[species_id %in% unlist(species_list)]
  
  return(data.table(
    pct_range_pac = mean(subset_data$pct_range_pac, na.rm = TRUE),
    pct_range_atl = mean(subset_data$pct_range_atl, na.rm = TRUE),
    pct_range_ind = mean(subset_data$pct_range_ind, na.rm = TRUE)
  ))
}
pixel_ocean_pcts <- pixel_species[, calculate_ocean_pcts(species_list, species_ocean_wide), by = .(x, y)]
pixel_ocean_pcts[, likely_ocean := fifelse(
  pct_range_pac > pct_range_atl & pct_range_pac > pct_range_ind, "Pacific",
  fifelse(pct_range_atl > pct_range_pac & pct_range_atl > pct_range_ind, "Atlantic",
          fifelse(pct_range_ind > pct_range_pac & pct_range_ind > pct_range_atl, "Indian", 
                  "Undetermined"))
)]
pixel_summary <- merge(pixel_ocean_pcts, species_per_pixel, by = c("x", "y"))
comparison <- merge(pixel_summary, ocean_lookup, by = c("x", "y"), all.x = TRUE)
comparison[, needs_reassignment := FALSE] 
comparison[current_ocean == "Atlantic" & (likely_ocean == "Pacific" | likely_ocean == "Indian"), 
           needs_reassignment := TRUE]
comparison[(current_ocean == "Pacific" | current_ocean == "Indian") & likely_ocean == "Atlantic",
           needs_reassignment := TRUE]
#2 pixels identified as needing reassignment
pixels_to_update <- data.table(
  x = c(20, 33),
  y = c(-35.5, 30.5),
  new_ocean = c("Pacific", "Indian")) #currently both assigned as atlantic
for (i in 1:nrow(pixels_to_update)) {
  x_val <- pixels_to_update[i, x]
  y_val <- pixels_to_update[i, y]
  new_ocean <- pixels_to_update[i, new_ocean]
  
  allcoral_coord[x == x_val & y == y_val, ocean := new_ocean]
}


# more cleaning of species only present in indo-pacific, not atlantic
pacific_species <- unique(species_coord_ocean[ocean == "Pacific", .(species_id)])
atlantic_species <- unique(species_coord_ocean[ocean == "Atlantic", .(species_id)])
indian_species <- unique(species_coord_ocean[ocean == "Indian", .(species_id)])
pacific_atlantic <- merge(pacific_species, atlantic_species, by = "species_id")
pacific_indian <- merge(pacific_species, indian_species, by = "species_id")
atlantic_indian <- merge(atlantic_species, indian_species, by = "species_id")
all_oceans <- merge(pacific_atlantic, indian_species, by = "species_id")

#fix points that should be assigned pacific, not atlantic
species_to_fix <- all_oceans$species_id

pixel_ocean_lookup <- unique(allcoral_coord[, .(x, y, ocean)])
setkey(pixel_ocean_lookup, x, y)

changes_made <- list()

for (sp_id in species_to_fix) {
  # Get all pixels for this SINGLE species
  sp_pixels <- species_coord[species_id == sp_id, .(x, y)]
  
  # Add current ocean assignment
  sp_pixels <- merge(sp_pixels, pixel_ocean_lookup, by = c("x", "y"), all.x = TRUE)
  
  # Identify Atlantic pixels that need to be reassigned
  atlantic_pixels <- sp_pixels[ocean == "Atlantic"]
  
  if (nrow(atlantic_pixels) == 0) {
    cat("No Atlantic pixels found for species", sp_id, "\n")
    next
  }
  
  # Find non-Atlantic pixels for this species
  valid_pixels <- sp_pixels[ocean %in% c("Pacific", "Indian")]
  
  if (nrow(valid_pixels) == 0) {
    cat("Warning: No Pacific/Indian pixels found for species", sp_id, "- cannot determine nearest ocean\n")
    next
  }
  
  # For each Atlantic pixel, find the nearest valid pixel
  for (i in 1:nrow(atlantic_pixels)) {
    current_pixel <- atlantic_pixels[i]
    
    # Calculate distances to all valid pixels
    coords_query <- matrix(c(current_pixel$x, current_pixel$y), ncol = 2)
    coords_reference <- as.matrix(valid_pixels[, .(x, y)])
    
    # Handle cases where longitude wraps around
    distances <- rep(Inf, nrow(coords_reference))
    for (j in 1:nrow(coords_reference)) {
      # Calculate direct distance
      dx <- coords_reference[j, 1] - coords_query[1]
      dy <- coords_reference[j, 2] - coords_query[2]
      
      # Check if wrapping around in x direction would make it closer
      # Assuming 360 degree wrap-around
      if (abs(dx) > 180) {
        if (dx > 0) dx = dx - 360
        else dx = dx + 360
      }
      
      distances[j] <- sqrt(dx^2 + dy^2)
    }
    
    # Find the index of the closest valid pixel
    closest_idx <- which.min(distances)
    new_ocean <- valid_pixels[closest_idx, ocean]
    
    # Record the change
    changes_made[[length(changes_made) + 1]] <- data.table(
      species_id = sp_id,
      x = current_pixel$x,
      y = current_pixel$y,
      old_ocean = "Atlantic",
      new_ocean = new_ocean
    )
    
    # Update the pixel ocean assignment in allcoral_coord
    allcoral_coord[x == current_pixel$x & y == current_pixel$y, ocean := new_ocean]
    
    # Update ocean_envelopes if it has the same structure
    if ("x" %in% names(ocean_envelopes) && "y" %in% names(ocean_envelopes) && "ocean" %in% names(ocean_envelopes)) {
      ocean_envelopes[x == current_pixel$x & y == current_pixel$y, ocean := new_ocean]
    }
  }
  
  cat("Fixed", nrow(atlantic_pixels), "Atlantic pixels for species", sp_id, "\n")
}

all_changes <- rbindlist(changes_made)


#fix points that should be assigned atlantic, not pacific
species_to_fix <- pacific_atlantic$species_id

pixel_ocean_lookup <- unique(allcoral_coord[, .(x, y, ocean)])
setkey(pixel_ocean_lookup, x, y)

problematic_pixels <- data.table()
for (sp_id in species_to_fix) {
  # Get all pixels for this species
  sp_pixels <- species_coord[species_id == sp_id, .(x, y)]
  
  # Add current ocean assignment
  sp_pixels <- merge(sp_pixels, pixel_ocean_lookup, by = c("x", "y"), all.x = TRUE)
  
  # Identify Pacific pixels that need to be reassigned to Atlantic
  pacific_pixels <- sp_pixels[ocean == "Pacific", .(species_id = sp_id, x, y, ocean)]
  
  if (nrow(pacific_pixels) > 0) {
    problematic_pixels <- rbind(problematic_pixels, pacific_pixels)
    cat("Found", nrow(pacific_pixels), "Pacific pixels for Atlantic-only species", sp_id, "\n")
  }
}
unique_xy <- unique(problematic_pixels[, c("x", "y")]) #three panama-adjacent pixels being assigned as both atlantic & pacific

problem_pxl_mask <- allcoral_coord$x %in% c(278, 279,281, 282) & 
  allcoral_coord$y %in% c(8.5, 9.5)

allcoral_coord$ocean[problem_pxl_mask] <- NA




#calculate envelopes
oceans <- unique(allcoral_coord$ocean)
ocean_envelopes <- data.table()

for (ocean_name in oceans) {
  if (is.na(ocean_name)) next
  
  ocean_data <- allcoral_coord[ocean == ocean_name]
  if (nrow(ocean_data) > 0) {
    ocean_envelope <- calculate_envelope(ocean_data, "allcoral", ocean_name)
    ocean_envelopes <- rbind(ocean_envelopes, ocean_envelope)
  }
}

setnames(ocean_envelopes, "realm", "ocean")
fwrite(ocean_envelopes, paste0(output_directory, 'ocean_envelopes_1982-1992.csv'))
fwrite(allcoral_coord, paste0(output_directory, 'allcoral_pixels_with_oceans_1982-1992.csv'))


#calculate envelope for indo-pacific
indo_pacific_data <- allcoral_coord[ocean %in% c("Indian", "Pacific")]
indo_pacific_envelope <- calculate_envelope(indo_pacific_data, "allcoral", "indo-pacific")
setnames(indo_pacific_envelope, "realm", "ocean")
ocean_envelopes <- rbind(ocean_envelopes, indo_pacific_envelope)


fwrite(ocean_envelopes, paste0(output_directory, 'ocean_envelopes_1982-1992.csv'))










# Add aragonite analysis



arag <- rast(paste0(cmip_folder, "arag/Aragonite_median_historical_foc.nc"))

# Extract layer 16 which corresponds to 1982-1992 period
arag_8292 <- arag[[16]]

# Rotate to match the coordinate system of other variables
arag180 <- rotate(arag_8292)
arag180r <- crop(arag180, ext(-180.5, -180, -90, 90))
arag180l <- crop(arag180, ext(-180, 179.5, -90, 90))
arag180r <- terra::shift(arag180r, dx = 360)
arag180 <- merge(arag180l, arag180r)

# Extract aragonite values for species coordinates
coords <- unique(species_coord[, .(x, y)])
points <- vect(coords, geom=c("x", "y"), crs=crs(arag180))

# Extract aragonite values at each point
arag_values <- terra::extract(arag180, points)

# Create a lookup table with coordinates and aragonite values
arag_lookup <- data.table(
  x = coords$x,
  y = coords$y,
  aragonite = arag_values[[2]])

# Merge aragonite values into the species coordinates data
allcoral_coord[, aragonite := NULL]
species_coord[, aragonite := NULL]
species_coord <- merge(species_coord, arag_lookup, by=c("x", "y"), all.x=TRUE)
#allcoral_coord <- merge(allcoral_coord, arag_lookup, by=c("x", "y"), all.x=TRUE)
allcoral_coord_temp <- copy(allcoral_coord)
allcoral_coord_temp[, x_adj := ifelse(x > 180, x - 360, x)]
allcoral_coord_temp <- merge(allcoral_coord_temp, 
                             arag_lookup[, .(x, y, aragonite)], 
                             by.x=c("x_adj", "y"), by.y=c("x", "y"), all.x=TRUE)
allcoral_coord[, aragonite := allcoral_coord_temp$aragonite]

# Calculate species-level aragonite statistics
species_arag <- species_coord[, .(
  aragonite_mean = mean(aragonite, na.rm=TRUE),
  aragonite_min = quantile(aragonite, probs = 0.025, na.rm = TRUE),
  aragonite_max = quantile(aragonite, probs = 0.975, na.rm = TRUE),
  aragonite_min_absolute = min(aragonite, na.rm = TRUE),
  aragonite_max_absolute = max(aragonite, na.rm = TRUE),
  aragonite_median = median(aragonite, na.rm = TRUE)
), by=species_id]

allcoral_arag <- allcoral_coord[, .(
  aragonite_mean = mean(aragonite, na.rm=TRUE),
  aragonite_min = quantile(aragonite, probs = 0.025, na.rm = TRUE),
  aragonite_max = quantile(aragonite, probs = 0.975, na.rm = TRUE),
  aragonite_min_absolute = min(aragonite, na.rm = TRUE),
  aragonite_max_absolute = max(aragonite, na.rm = TRUE),
  aragonite_median = median(aragonite, na.rm = TRUE)
), by=species_id]

# Merge aragonite statistics into the species envelope data
species_arag$species_id <- as.character(species_arag$species_id)
allcoral_arag$species_id <- as.character(allcoral_arag$species_id)
species_envelope$id_no <- as.character(species_envelope$id_no)
species_envelope <- merge(species_envelope, species_arag, by.x="id_no", by.y="species_id", all.x=TRUE)
allcoral_envelope <- merge(allcoral_envelope, allcoral_arag, by.x="id_no", by.y="species_id", all.x=TRUE)

# Save updated data with aragonite information
fwrite(species_coord, paste0(output_directory, 'species_pixels_1982-1992_with_arag.csv'))
fwrite(species_envelope, paste0(output_directory, 'climate_envelopes_1982-1992_with_arag.csv'))
fwrite(allcoral_envelope, paste0(output_directory, 'allcoral_envelopes_1982-1992_with_arag.csv'))
fwrite(allcoral_coord, paste0(output_directory, 'allcoral_pixels_1982-1992_with_arag.csv'))

