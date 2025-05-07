# figures
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)


source("functions.R")

scenarios <- c("ssp126","ssp245","ssp370","ssp585")

#open envelope exceedance data
exceedsp126 <- fread(paste0(output_directory, "species_exceedance_summary_", scenarios[[1]], ".csv"))
exceedsp245 <- fread(paste0(output_directory, "species_exceedance_summary_", scenarios[[2]], ".csv"))
exceedsp370 <- fread(paste0(output_directory, "species_exceedance_summary_", scenarios[[3]], ".csv"))
exceedsp585 <- fread(paste0(output_directory, "species_exceedance_summary_", scenarios[[4]], ".csv"))
exceedall126 <- fread(paste0(output_directory, "allcoral_exceedance_summary_", scenarios[[1]], ".csv"))
exceedall245 <- fread(paste0(output_directory, "allcoral_exceedance_summary_", scenarios[[2]], ".csv"))
exceedall370 <- fread(paste0(output_directory, "allcoral_exceedance_summary_", scenarios[[3]], ".csv"))
exceedall585 <- fread(paste0(output_directory, "allcoral_exceedance_summary_", scenarios[[4]], ".csv"))
#open dhw data
dhwsp126 <- fread(paste0(output_directory, "coralsp_dhwsummary_annual_", scenarios[[1]], ".csv"))
dhwsp245 <- fread(paste0(output_directory, "coralsp_dhwsummary_annual_", scenarios[[2]], ".csv"))
dhwsp370 <- fread(paste0(output_directory, "coralsp_dhwsummary_annual_", scenarios[[3]], ".csv"))
dhwsp585 <- fread(paste0(output_directory, "coralsp_dhwsummary_annual_", scenarios[[4]], ".csv"))
dhwall126 <- fread(paste0(output_directory, "allcoral_dhwsummary_annual_", scenarios[[1]], ".csv"))
dhwall245 <- fread(paste0(output_directory, "allcoral_dhwsummary_annual_", scenarios[[2]], ".csv"))
dhwall370 <- fread(paste0(output_directory, "allcoral_dhwsummary_annual_", scenarios[[3]], ".csv"))
dhwall585 <- fread(paste0(output_directory, "allcoral_dhwsummary_annual_", scenarios[[4]], ".csv"))
dhwall585$total_area <- exceedall585$total_area

allcoral_envelope <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992.csv'))
species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))

exceed585 <- rbind(exceedsp585,exceedall585)
exceedall585 <- merge(exceedall585,allcoral_envelope,by.x="species_id",by.y="id_no",all.x=T)
species_envelope$id_no <- as.character(species_envelope$id_no)
exceedsp585$species_id <- as.character(exceedsp585$species_id)
exceedsp585 <- merge(exceedsp585,species_envelope,by.x="species_id",by.y="id_no",all.x=T)
dhwsp585$species_id <- as.character(dhwsp585$species_id)
dhwsp585 <- merge(dhwsp585,species_envelope,by.x="species_id",by.y="id_no",all.x=T)



### FIG 1
#temperature envelope exceeded, bleaching events, and length of bleaching. pH envelope exceeded

exceed585$pct_exceed_temp_area <- exceed585$temp_exceed_area / exceed585$total_area
exceed585$pct_exceed_ph_area <- exceed585$ph_exceed_area / exceed585$total_area
exceed585$pct_exceed_dual_area <- exceed585$dual_exceed_area / exceed585$total_area
exceed585$pct_dual_or_exceed_area <- exceed585$dual_or_exceed_area / exceed585$total_area



### A

ocean_data <- fread(paste0(output_directory, "all_oceans_exceedance_summary_ssp585.csv"))
ocean_data[, pct_exceed_temp_area := temp_exceed_area / total_area]


# percentage area of overall coral exceeded envelope vs percent area of species exceeded envelopes
temp_envelope <- ggplot() +
  # average of all species 
  geom_line(data = exceed585_all_oceans %>% 
              filter(species_id != "allcoral", ocean == "indo-pacific") %>%
              group_by(year, scenario) %>%
              summarize(mean_unsuitable = mean(pct_exceed_temp_area*100, na.rm = TRUE), .groups = "drop"),
            aes(x = year, y = mean_unsuitable, group = scenario, 
                color = "Indo-Pacific"),linewidth = 1, linetype = "dashed") +
  geom_line(data = exceed585_all_oceans %>% 
              filter(species_id != "allcoral", ocean == "Atlantic") %>%
              group_by(year, scenario) %>%
              summarize(mean_unsuitable = mean(pct_exceed_temp_area*100, na.rm = TRUE), .groups = "drop"),
            aes(x = year, y = mean_unsuitable, group = scenario,
                color = "Atlantic"),linewidth = 1, linetype = "dashed") +
  # Add ocean-specific lines
  geom_line(data = ocean_data %>% filter(ocean == "indo-pacific"),
            aes(x = year, y = pct_exceed_temp_area*100, group = ocean,
                color = "Indo-Pacific"),linewidth = 1, linetype = "solid") +
  geom_line(data = ocean_data %>% filter(ocean == "Atlantic"),
            aes(x = year, y = pct_exceed_temp_area*100, group = ocean,
                color = "Atlantic"),linewidth = 1, linetype = "solid") +
  scale_color_manual(name = "Ocean",
                     values = c("Atlantic" = "blue", "Indo-Pacific" = "red")) +
  scale_linetype_manual(name = "Type",
                        values = c("Species mean" = "dashed", "Ecosystem mean" = "solid")) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                     breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Year",
       y = "Habitat loss (area exceeding \ntemperature envelopes (%))") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = c(0.95, 0.05),  # x,y coordinates (0.5 is center, 0.1 is near bottom)
    legend.justification = c(0.5, 0), # anchor point for the legend
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5),
    legend.margin = margin(6, 6, 6, 6))



### B
# mean temperature increase in ranges

exceedall585$mean_temp <- exceedall585$mean_temp_diff + exceedall585$sst_max
exceedall585$mean_temp_delta <- exceedall585$mean_temp - exceedall585$sst_summer_mean 

exceedsp585$mean_temp <- exceedsp585$mean_temp_diff + exceedsp585$sst_max
exceedsp585$mean_temp_delta <- exceedsp585$mean_temp - exceedsp585$sst_summer_mean 



temp_delta <-  ggplot() +
  #species lines
  geom_line(data = exceedsp585 %>% filter(species_id != "allcoral"),
            aes(x = year, 
                y = mean_temp_delta,
                group = interaction(species_id, scenario)), 
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceedsp585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_summer_delta = mean(mean_temp_delta, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_summer_delta,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceedall585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = mean_temp_delta,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(-0.1, 6),breaks=(seq(0, 6, 1)),expand=c(0,0)) +
  scale_x_continuous(limits = c(2015, 2100),breaks=(seq(2020, 2100, 20)),expand=c(0,0))+
  labs(x = "Year",
       y = "Exposure (Summer temperature \ndifference (°C))") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )





#temperature and pH profiles

calculate_climate_profile <- function(output_directory, scenario, cmip_folder) {
  # Read current coral data with summer conditions
  allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))
  
  # Calculate current summer means (1982-1992)
  # Extract summer months from your existing data (sst_8292 and ph_8292)
  sst <- rast(paste0(cmip_folder,'tos/tos_Omon_modelmeanfoc_hist_185001-201412.tif'))
  ph <- rast(paste0(cmip_folder,'ph/ph_Omon_modelmeanfoc_hist_185001-201412.tif'))
  
  # Get the indices for 1982-1992
  start_8292 <- (1982-1850)*12 + 1
  end_8292 <- (1992-1850)*12 + 12
  sst_8292 <- sst[[start_8292:end_8292]]
  ph_8292 <- ph[[start_8292:end_8292]]
  
  # Extract summer months
  # Northern hemisphere summer (July-September)
  years <- 11  # 1982-1992
  julyaugsept <- c(7, 8, 9)
  julyaugsepts <- numeric(0)
  for (i in 0:(years - 1)) {
    julyaugsepts <- c(julyaugsepts, julyaugsept + i * 12)
  }
  sst_summerN_8292 <- subset(sst_8292, julyaugsepts)
  ph_summerN_8292 <- subset(ph_8292, julyaugsepts)
  
  # Southern hemisphere summer (January-March)
  janfebmarch <- c(1, 2, 3)
  janfebmarchs <- numeric(0)
  for (i in 0:(years - 1)) {
    janfebmarchs <- c(janfebmarchs, janfebmarch + i * 12)
  }
  sst_summerS_8292 <- subset(sst_8292, janfebmarchs)
  ph_summerS_8292 <- subset(ph_8292, janfebmarchs)
  
  # Calculate mean for summer months
  sst_summerN_mean_8292 <- app(sst_summerN_8292, fun=mean, na.rm=TRUE)
  sst_summerS_mean_8292 <- app(sst_summerS_8292, fun=mean, na.rm=TRUE)
  ph_summerN_mean_8292 <- app(ph_summerN_8292, fun=mean, na.rm=TRUE)
  ph_summerS_mean_8292 <- app(ph_summerS_8292, fun=mean, na.rm=TRUE)
  
  # Read future summer rasters
  tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OyrsummerN_modelmeanfoc_',scenario,'_201501-210012.tif'))
  ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OyrsummerS_modelmeanfoc_',scenario,'_201501-210012.tif'))
  
  # Extract the 86th layer (year 2100) from the future rasters
  tos_summerN_2100 <- tos_summerN[[86]]
  tos_summerS_2100 <- tos_summerS[[86]]
  ph_summerN_2100 <- ph_summerN[[86]]
  ph_summerS_2100 <- ph_summerS[[86]]
  
  # Extract values at coordinates
  # Determine hemisphere for each coordinate
  allcoral_coord[, hemisphere := ifelse(y > 0, "N", "S")]
  
  # Extract current SST values based on hemisphere
  north_coords <- allcoral_coord[hemisphere == "N"]
  south_coords <- allcoral_coord[hemisphere == "S"]
  
  if (nrow(north_coords) > 0) {
    sst_north_values <- terra::extract(sst_summerN_mean_8292, north_coords[, .(x, y)])
    ph_north_values <- terra::extract(ph_summerN_mean_8292, north_coords[, .(x, y)])
    allcoral_coord[hemisphere == "N", sst_summer := sst_north_values[,2]]
    allcoral_coord[hemisphere == "N", ph_summer := ph_north_values[,2]]
  }
  
  if (nrow(south_coords) > 0) {
    sst_south_values <- terra::extract(sst_summerS_mean_8292, south_coords[, .(x, y)])
    ph_south_values <- terra::extract(ph_summerS_mean_8292, south_coords[, .(x, y)])
    allcoral_coord[hemisphere == "S", sst_summer := sst_south_values[,2]]
    allcoral_coord[hemisphere == "S", ph_summer := ph_south_values[,2]]
  }
  
  # Extract future SST and pH values based on hemisphere
  if (nrow(north_coords) > 0) {
    sst_north_future <- terra::extract(tos_summerN_2100, north_coords[, .(x, y)])
    ph_north_future <- terra::extract(ph_summerN_2100, north_coords[, .(x, y)])
    allcoral_coord[hemisphere == "N", sst_future := sst_north_future[,2]]
    allcoral_coord[hemisphere == "N", ph_future := ph_north_future[,2]]
  }
  
  if (nrow(south_coords) > 0) {
    sst_south_future <- terra::extract(tos_summerS_2100, south_coords[, .(x, y)])
    ph_south_future <- terra::extract(ph_summerS_2100, south_coords[, .(x, y)])
    allcoral_coord[hemisphere == "S", sst_future := sst_south_future[,2]]
    allcoral_coord[hemisphere == "S", ph_future := ph_south_future[,2]]
  }
  
  # Calculate current latitude data
  current_lat <- allcoral_coord[, .(
    temperature = mean(sst_summer, na.rm = TRUE),
    ph = mean(ph_summer, na.rm = TRUE)
  ), by = y]
  setnames(current_lat, "y", "latitude")
  current_lat[, scenario := "Today"]
  
  # Calculate future latitude data
  future_lat <- allcoral_coord[, .(
    temperature = mean(sst_future, na.rm = TRUE),
    ph = mean(ph_future, na.rm = TRUE)
  ), by = y]
  setnames(future_lat, "y", "latitude")
  future_lat[, scenario := "Future"]
  
  # Combine current and future data
  combined_lat <- rbind(current_lat, future_lat)
  
  return(combined_lat)
}

for(scenario in scenarios){
  var_name <- paste0('combined_lat',gsub('ssp','',scenario))
  assign(var_name, calculate_climate_profile(output_directory, scenario,cmip_folder), envir = .GlobalEnv)
}

#temperature profile
temp_lat <- ggplot(combined_lat585, aes(x = latitude, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red"),
                     labels = c("Today" = "Today", "Future" = "SSP585")) +
  theme_classic() +
  labs(x = "Latitude (degrees)", 
       y = "Mean summer \ntemperature (°C)",
       color = "Scenario") +
  scale_x_continuous(limits=c(-45,45),expand=c(0,0),breaks = seq(-45, 45, 45)) +
  scale_y_continuous(limits=c(15,35),breaks = seq(15, 35, 5),expand=c(0,0))+
  theme(legend.position = c(0.5, 0.1),  # x,y coordinates (0.5 is center, 0.1 is near bottom)
        legend.justification = c(0.5, 0), # anchor point for the legend
        legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5),
        legend.margin = margin(6, 6, 6, 6))



#temp increase v range size
tempdelta_range <- ggplot() +
  geom_pointdensity(data = exceedsp585[exceedsp585$year==2100,], 
                    aes(x = total_area/exceedall585$total_area[1]*100, y = mean_temp_delta)) +
  scale_color_gradient(low='#cbecff', high="#2171b5") +
  geom_smooth(data = exceedsp585[exceedsp585$year==2100,],
              aes(x = total_area/exceedall585$total_area[1]*100, y = mean_temp_delta),
              method = "loess", span = 0.3, color = "black", se=F) +
  geom_point(data = exceedall585[exceedall585$year==2100,], 
             aes(x = 100, y = mean_temp_delta), 
             color = "red", 
             size = 3) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain area (%)", 
       y = "Temperature exposure \nby 2100 (°C)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,3))+
  scale_y_continuous(expand=c(0,0),
                     limits=c(3,6),
                     breaks=seq(3,6,0.5))+
  theme(legend.position="none")


# % loss vs range size
# Sort the data by range size
sorted_data <- exceedsp585[exceedsp585$year==2100,] %>%
  arrange(total_area/exceedall585$total_area[1])
# Calculate cumulative percentage of species
sorted_data$cum_species_pct <- seq_len(nrow(sorted_data)) / nrow(sorted_data)

#plot3 with critical range size shading
displacement <- calculate_latitudinal_displacement(combined_lat585 %>% filter(scenario == "Today") %>% select(latitude,temperature),
                                                   combined_lat585 %>% filter(scenario == "Future") %>% select(latitude,temperature))
# Filter out cross-hemisphere matches
displacement <- displacement %>%
  mutate(cross_hemisphere = !is.na(latitude) & !is.na(future_lat) & sign(latitude) != sign(future_lat)) %>%
  # Set displacement to NA for cross-hemisphere matches
  mutate(displacement = ifelse(cross_hemisphere, NA, displacement),
    display_displacement = ifelse(cross_hemisphere, NA, display_displacement))

min_critical_range <- min(displacement$displacement, na.rm = TRUE)
max_critical_range <- max(displacement$displacement, na.rm = TRUE)

tempexceed_range <- ggplot(exceed585[exceed585$year==2100,] %>% filter(species_id != "allcoral"), 
                           aes(x = total_area/exceedall585$total_area[1]*100, y = pct_exceed_temp_area*100)) +
  geom_pointdensity()+
  scale_color_gradient(low='#cbecff',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  # Add the cumulative species line
  geom_line(data = sorted_data, aes(x = total_area/exceedall585$total_area[1]*100, y = cum_species_pct*100), color = "blue", linewidth = 1) +
  labs(x = "Range as proportion \nof domain size", 
       y = "Area habitat loss \nby 2100 (%)")+
  scale_y_continuous(
    name = "Area exceeding\n thermal envelope (%)",
    sec.axis = sec_axis(~. * 100 / max(results_sim$exceed_max_pct), 
                        name = "Cumulative species (%)"),
    expand=c(0,0),
    limits = c(0, 100.5),
    breaks = seq(0, 100, 20)
  ) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain area (%)", 
       y = "Area exceeding envelope (%)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  theme(legend.position="none",
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.line.y.right = element_line(color="blue"))
  








# % of species with range loss greater than the largest-range species loss in 2100
exceed585_2100 <- exceed585[exceed585$year == 2100, ]
exceed585_2100 <- exceed585_2100 %>% filter(species_id != "allcoral")
exceed585_2100$pct_range_domain <- exceed585_2100$total_area / exceedall585$total_area[1] * 100
bins <- seq(0, 100, by = 5)
exceed585_2100$bin <- cut(exceed585_2100$pct_range_domain, 
                       breaks = bins,
                       labels = paste0(bins[-length(bins)], "-", bins[-1], "%"),
                       include.lowest = TRUE,
                       right = FALSE)
largest_bin_data <- exceed585_2100[exceed585_2100$bin == "65-70%", ]
largest_bin_mean_loss <- mean(largest_bin_data$pct_exceed_temp_area, na.rm = TRUE)
bin_results <- data.frame()
for (b in levels(exceed585_2100$bin)) {
  bin_data <- exceed585_2100[exceed585_2100$bin == b, ]
  total_species <- nrow(bin_data)
  
  if (total_species > 0) {
    # Count species exceeding the threshold
    exceeding_threshold <- sum(bin_data$pct_exceed_temp_area > largest_bin_mean_loss)
    percentage <- (exceeding_threshold / total_species) * 100
    
    # Calculate standard error (for error bars)
    p <- percentage / 100
    se <- sqrt((p * (1 - p)) / total_species) * 100
    
    # Add to results
    bin_results <- rbind(bin_results, data.frame(
      bin = b,
      bin_center = mean(as.numeric(gsub("%", "", unlist(strsplit(gsub("-", " ", b), " "))))),
      percentage = percentage,
      se = se,
      count = total_species
    ))
  }
}

pct_loss_compare_largest <- ggplot(bin_results, aes(x = bin_center, y = percentage)) +
  # Regular points
  geom_point(data = bin_results[bin_results$bin_center != 67.5, ], size = 3) +
  # Reference point in red
  geom_point(data = bin_results[bin_results$bin_center == 67.5, ], size = 3, color = "red") +
  geom_errorbar(aes(ymin = pmax(0, percentage - se), 
                    ymax = pmin(100, percentage + se)), 
                width = 1.5) +
  # More concise y-axis title
  labs(x = "Range size (% of domain area)",
       y = "% species exceeding \nlargest-range habitat loss") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, by = 20),
                     limits = c(0, 105),
                     expand = c(0, 0)) +
  # Optional: add annotation explaining the reference point
  annotate("text", x = 83, y = 60, 
           label = "Reference \nthreshold", 
           color = "red", size = 3.5)














#distribution plot - sizes of ranges
range_size <- ggplot(exceedsp585[exceedsp585$year==2100,],aes(x=total_area/exceedall585$total_area[1]*100))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Range as proportion \nof domain area (%)", 
       y = "Number of species")+
  scale_x_continuous(limits = c(-3, 100), 
                     breaks = seq(0, 100, 20),expand=c(0,0))+
  scale_y_continuous(limits=c(0,150),expand=c(0,0))









### 
# range size v distance to analog

species_coord <- fread(paste0(output_directory, 'species_pixels_1982-1992.csv'))

# Bounding box approach for calculating species range in kilometers
species_range <- species_coord %>%
  group_by(species_id) %>%
  summarize(
    # Store original min/max for reference
    x_min = min(x),
    x_max = max(x),
    
    # Latitudinal range (y coordinates)
    y_min = min(y),
    y_max = max(y),
    y_range = max(y) - min(y),
    
    # For longitude range calculation, handle the dateline crossing
    .x_360 = list(ifelse(x < 0, x + 360, x))
  ) %>%
  rowwise() %>%
  mutate(
    # Calculate longitude range accounting for dateline
    x_range_corrected = {
      # Sort the longitudes
      x_sorted <- sort(unlist(.x_360))
      
      # Find the maximum gap between consecutive points
      n <- length(x_sorted)
      if (n <= 1) {
        return(0)
      }
      
      gaps <- c(diff(x_sorted), (x_sorted[1] + 360) - x_sorted[n])
      max_gap <- max(gaps)
      
      # The actual range is 360 minus the maximum gap
      360 - max_gap
    },
    
    # Add the maximum of latitudinal and longitudinal extents
    max_latlon_range = max(y_range, x_range_corrected),
    
    # Calculate range in kilometers using a bounding box approach
    range_km = {
      # Create an sf polygon representing the bounding box
      # Need to handle dateline crossing
      if (x_range_corrected > 180) {
        # For species that cross the dateline, we need a different approach
        # Convert all coordinates to the 0-360 system
        x_360 <- unlist(.x_360)
        min_x_360 <- min(x_360)
        max_x_360 <- max(x_360)
        
        # Convert back to -180 to 180 system
        min_x_180 <- ifelse(min_x_360 > 180, min_x_360 - 360, min_x_360)
        max_x_180 <- ifelse(max_x_360 > 180, max_x_360 - 360, max_x_360)
        
        # Create bounding box
        bbox <- st_bbox(c(xmin = min_x_180, ymin = y_min, 
                          xmax = max_x_180, ymax = y_max), 
                        crs = 4326)
      } else {
        # Standard case - no dateline crossing
        bbox <- st_bbox(c(xmin = x_min, ymin = y_min, 
                          xmax = x_max, ymax = y_max), 
                        crs = 4326)
      }
      
      # Convert to polygon
      bbox_poly <- st_as_sfc(bbox)
      
      # Calculate maximum distance across the bounding box (diagonal)
      pts <- st_cast(bbox_poly, "POINT")
      dist_matrix <- st_distance(pts)
      max(dist_matrix) / 1000  # Convert to kilometers
    },
    
    # Add calculations for latitudinal and longitudinal ranges in km
    lat_range_km = {
      # Create points at the same longitude but different latitudes
      p1 <- st_point(c(x_min, y_min))
      p2 <- st_point(c(x_min, y_max))
      
      # Convert to sf objects with CRS
      p1_sf <- st_sfc(p1, crs = 4326)
      p2_sf <- st_sfc(p2, crs = 4326)
      
      # Calculate distance
      as.numeric(st_distance(p1_sf, p2_sf)) / 1000  # Convert to kilometers
    },
    
    lon_range_km = {
      # For longitude, we need to measure at the same latitude
      # Use the average latitude for better accuracy
      avg_lat <- (y_min + y_max) / 2
      
      if (x_range_corrected > 180) {
        # Handle dateline crossing
        x_360 <- unlist(.x_360)
        min_x_360 <- min(x_360)
        max_x_360 <- max(x_360)
        
        # Convert back to -180 to 180 system
        min_x_180 <- ifelse(min_x_360 > 180, min_x_360 - 360, min_x_360)
        max_x_180 <- ifelse(max_x_360 > 180, max_x_360 - 360, max_x_360)
        
        p1 <- st_point(c(min_x_180, avg_lat))
        p2 <- st_point(c(max_x_180, avg_lat))
      } else {
        p1 <- st_point(c(x_min, avg_lat))
        p2 <- st_point(c(x_max, avg_lat))
      }
      
      # Convert to sf objects with CRS
      p1_sf <- st_sfc(p1, crs = 4326)
      p2_sf <- st_sfc(p2, crs = 4326)
      
      # Calculate distance
      as.numeric(st_distance(p1_sf, p2_sf)) / 1000  # Convert to kilometers
    }
  ) %>%
  ungroup() %>%  # Ungroup before using select
  dplyr::select(-starts_with("."))  # Use dplyr::select explicitly and remove temporary columns


diffanalogN_585 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerN', "_", "ssp585", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogS_585 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerS', "_", "ssp585", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogN_585 <- diffanalogN_585[, .SD[1], by = .(px_temp, py_temp)]
diffanalogS_585 <- diffanalogS_585[, .SD[1], by = .(px_temp, py_temp)]


species_coord$px_temp <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
species_coord[, py_temp := y] 

north_result <- merge(species_coord[y >= 0], diffanalogN_585, by = c("px_temp", "py_temp"), all.x = TRUE)
south_result <- merge(species_coord[y < 0], diffanalogS_585, by = c("px_temp", "py_temp"), all.x = TRUE)
species_analog585 <- rbindlist(list(north_result, south_result), fill = TRUE)

sp_rangesize <- unique(exceedsp585[, .(species_id, total_area)])
sp_rangesize$species_id <- as.character(sp_rangesize$species_id)
species_analog585$species_id <- as.character(species_analog585$species_id)
species_analog585 <- merge(species_analog585, sp_rangesize, by = "species_id", all.x = TRUE)

species_avg <- species_analog585[, .(
  avg_distance_temp = mean(distance_temp, na.rm = TRUE),
  se_distance_temp = sd(distance_temp, na.rm = TRUE) / sqrt(.N),
  avg_distance_ph = mean(distance_ph, na.rm = TRUE),
  se_distance_ph = sd(distance_ph, na.rm = TRUE) / sqrt(.N),
  total_area = first(total_area)  # Take the first total_area value for each species
), by = species_id]











### SSP distance to analog v % complete habitat loss

# calculate summary data for analog distances per SSP
analog_stats_585 <- data.frame(
  ssp = "ssp585",
  mean_dist = mean(species_avg$avg_distance_temp, na.rm = TRUE),
  min_dist = min(species_avg$avg_distance_temp, na.rm = TRUE),
  max_dist = max(species_avg$avg_distance_temp, na.rm = TRUE))

# Process SSP126 analog distances 
diffanalogN_126 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerN', "_", "ssp126", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogS_126 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerS', "_", "ssp126", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogN_126 <- diffanalogN_126[, .SD[1], by = .(px_temp, py_temp)]
diffanalogS_126 <- diffanalogS_126[, .SD[1], by = .(px_temp, py_temp)]

# Merge with species coordinates (similar to what you did for SSP585)
north_result_126 <- merge(species_coord[y >= 0], diffanalogN_126, by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
south_result_126 <- merge(species_coord[y < 0], diffanalogS_126,  by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
species_analog126 <- rbindlist(list(north_result_126, south_result_126), fill = TRUE)
species_analog126$species_id <- as.character(species_analog126$species_id)
species_analog126 <- merge(species_analog126, sp_rangesize, by = "species_id", all.x = TRUE)

# Calculate stats per species
species_avg_126 <- species_analog126[, .(
  avg_distance_temp = mean(distance_temp, na.rm = TRUE),
  se_distance_temp = sd(distance_temp, na.rm = TRUE) / sqrt(.N),
  avg_distance_ph = mean(distance_ph, na.rm = TRUE),
  se_distance_ph = sd(distance_ph, na.rm = TRUE) / sqrt(.N),
  total_area = first(total_area)
), by = species_id]

# Add to stats dataframe
analog_stats_126 <- data.frame(
  ssp = "ssp126",
  mean_dist = mean(species_avg_126$avg_distance_temp, na.rm = TRUE),
  min_dist = min(species_avg_126$avg_distance_temp, na.rm = TRUE),
  max_dist = max(species_avg_126$avg_distance_temp, na.rm = TRUE)
)

# Process SSP245 analog distances 
diffanalogN_245 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerN', "_", "ssp245", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogS_245 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerS', "_", "ssp245", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogN_245 <- diffanalogN_245[, .SD[1], by = .(px_temp, py_temp)]
diffanalogS_245 <- diffanalogS_245[, .SD[1], by = .(px_temp, py_temp)]

# Merge with species coordinates (similar to what you did for SSP585)
north_result_245 <- merge(species_coord[y >= 0], diffanalogN_245,  by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
south_result_245 <- merge(species_coord[y < 0], diffanalogS_245,  by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
species_analog245 <- rbindlist(list(north_result_245, south_result_245), fill = TRUE)
species_analog245$species_id <- as.character(species_analog245$species_id)
species_analog245 <- merge(species_analog245, sp_rangesize, by = "species_id", all.x = TRUE)

# Calculate stats per species
species_avg_245 <- species_analog245[, .(
  avg_distance_temp = mean(distance_temp, na.rm = TRUE),
  se_distance_temp = sd(distance_temp, na.rm = TRUE) / sqrt(.N),
  avg_distance_ph = mean(distance_ph, na.rm = TRUE),
  se_distance_ph = sd(distance_ph, na.rm = TRUE) / sqrt(.N),
  total_area = first(total_area)
), by = species_id]

# Add to stats dataframe
analog_stats_245 <- data.frame(
  ssp = "ssp245",
  mean_dist = mean(species_avg_245$avg_distance_temp, na.rm = TRUE),
  min_dist = min(species_avg_245$avg_distance_temp, na.rm = TRUE),
  max_dist = max(species_avg_245$avg_distance_temp, na.rm = TRUE)
)

# Process SSP370 analog distances 
diffanalogN_370 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerN', "_", "ssp370", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogS_370 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerS', "_", "ssp370", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogN_370 <- diffanalogN_370[, .SD[1], by = .(px_temp, py_temp)]
diffanalogS_370 <- diffanalogS_370[, .SD[1], by = .(px_temp, py_temp)]

# Merge with species coordinates (similar to what you did for SSP585)
north_result_370 <- merge(species_coord[y >= 0], diffanalogN_370,  by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
south_result_370 <- merge(species_coord[y < 0], diffanalogS_370,  by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
species_analog370 <- rbindlist(list(north_result_370, south_result_370), fill = TRUE)
species_analog370$species_id <- as.character(species_analog370$species_id)
species_analog370 <- merge(species_analog370, sp_rangesize, by = "species_id", all.x = TRUE)

# Calculate stats per species
species_avg_370 <- species_analog370[, .(
  avg_distance_temp = mean(distance_temp, na.rm = TRUE),
  se_distance_temp = sd(distance_temp, na.rm = TRUE) / sqrt(.N),
  avg_distance_ph = mean(distance_ph, na.rm = TRUE),
  se_distance_ph = sd(distance_ph, na.rm = TRUE) / sqrt(.N),
  total_area = first(total_area)
), by = species_id]

# Add to stats dataframe
analog_stats_370 <- data.frame(
  ssp = "ssp370",
  mean_dist = mean(species_avg_370$avg_distance_temp, na.rm = TRUE),
  min_dist = min(species_avg_370$avg_distance_temp, na.rm = TRUE),
  max_dist = max(species_avg_370$avg_distance_temp, na.rm = TRUE)
)

# Combine all analog stats
analog_stats <- rbind(analog_stats_126, analog_stats_245, analog_stats_370, analog_stats_585)

# Calculate percentage of species with complete habitat loss by 2100 for each SSP
# For SSP585
exceedsp585_2100 <- exceedsp585[year == 2100]
total_species_585 <- length(unique(exceedsp585_2100$species_id))
complete_loss_585 <- exceedsp585_2100[, .(
  habitat_loss_percent = (temp_exceed_area / total_area) * 100
), by = species_id]
complete_loss_count_585 <- nrow(complete_loss_585[habitat_loss_percent >= 99])
percent_loss_585 <- (complete_loss_count_585 / total_species_585) * 100

# SSP126
exceedsp126_2100 <- exceedsp126[year == 2100]
total_species_126 <- length(unique(exceedsp126_2100$species_id))
complete_loss_126 <- exceedsp126_2100[, .(
  habitat_loss_percent = (temp_exceed_area / total_area) * 100
), by = species_id]
complete_loss_count_126 <- nrow(complete_loss_126[habitat_loss_percent >= 99])
percent_loss_126 <- (complete_loss_count_126 / total_species_126) * 100

# SSP245 
exceedsp245_2100 <- exceedsp245[year == 2100]
total_species_245 <- length(unique(exceedsp245_2100$species_id))
complete_loss_245 <- exceedsp245_2100[, .(
  habitat_loss_percent = (temp_exceed_area / total_area) * 100
), by = species_id]
complete_loss_count_245 <- nrow(complete_loss_245[habitat_loss_percent >= 99])
percent_loss_245 <- (complete_loss_count_245 / total_species_245) * 100

# SSP370
exceedsp370_2100 <- exceedsp370[year == 2100]
total_species_370 <- length(unique(exceedsp370_2100$species_id))
complete_loss_370 <- exceedsp370_2100[, .(
  habitat_loss_percent = (temp_exceed_area / total_area) * 100
), by = species_id]
complete_loss_count_370 <- nrow(complete_loss_370[habitat_loss_percent >= 99])
percent_loss_370 <- (complete_loss_count_370 / total_species_370) * 100

# Combine habitat loss percentages
habitat_loss_data <- data.frame(
  ssp = c("ssp126", "ssp245", "ssp370", "ssp585"),
  percent_loss = c(percent_loss_126, percent_loss_245, percent_loss_370, percent_loss_585))

# Merge analog stats and habitat loss data
analog_habitatloss_data <- merge(analog_stats, habitat_loss_data, by = "ssp")

# Create more readable SSP labels
analog_habitatloss_data$ssp_label <- factor(analog_habitatloss_data$ssp, 
                              levels = c("ssp126", "ssp245", "ssp370", "ssp585"),
                              labels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))

# Order SSPs by warming scenario
analog_habitatloss_data <- analog_habitatloss_data[order(analog_habitatloss_data$ssp_label), ]

# Create the plot
analog_loss <- ggplot(analog_habitatloss_data, aes(x = mean_dist, y = percent_loss)) +
  geom_point(size = 3,alpha=0.5) +
  geom_segment(aes(x = min_dist, xend = max_dist, y = percent_loss, yend = percent_loss), 
               size = 1,alpha=0.5) +
  geom_text(aes(x = max_dist, label = ssp_label), 
            hjust = -0.1, 
            vjust = 0.3,  
            size = 3.5) +  
  scale_x_continuous(limits=c(0,4000),label=scales::comma,expand=c(0,0))+
  scale_y_continuous(limits=c(0,40),expand=c(0,0))+
  labs(x = "Distance to temperature analog (km)",
       y = "Species with complete \nhabitat loss by 2100 (%)") +
  theme_classic() 







species_range$species_id <- as.character(species_range$species_id)
species_range$range_km <- as.numeric(species_range$range_km)

species_avg <- merge(species_avg, exceedsp585 %>% dplyr::select(species_id, total_area ), by = "species_id", all.x = TRUE)
species_avg_370 <- merge(species_avg_370, exceedsp585 %>% dplyr::select(species_id, total_area ), by = "species_id", all.x = TRUE)
species_avg_245 <- merge(species_avg_245, exceedsp585 %>% dplyr::select(species_id, total_area ), by = "species_id", all.x = TRUE)
species_avg_126 <- merge(species_avg_126, exceedsp585 %>% dplyr::select(species_id, total_area ), by = "species_id", all.x = TRUE)

species_avg_126$scenario <- "SSP126"
species_avg_245$scenario <- "SSP245" 
species_avg_370$scenario <- "SSP370"
species_avg$scenario <- "SSP585"
species_avg_all<- rbind(species_avg_126, species_avg_245, species_avg_370,species_avg)

#just ssp585
ggplot(species_avg, aes(x = lat_range_km, y = avg_distance_temp)) +
  geom_point(size = 0.1, alpha = 0.2) +  # Small points with high transparency
  geom_errorbar(aes(ymin = avg_distance_temp - se_distance_temp, 
                    ymax = avg_distance_temp + se_distance_temp),
                width = 0.1, alpha = 0.2) +  # Error bars with transparency
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(label=scales::comma,limits=c(0,10000),breaks=c(seq(0,10000,5000)),expand=c(0,0))+
  scale_y_continuous(label=scales::comma,limits = c(0,3200),breaks=(seq(0, 3000, 1000)),expand=c(0,0))+
  labs(x = "Range size (km)",
       y = "Distance to \ntemperature analog (km)") +
  theme_classic() 


#all ssp
ggplot(species_avg_all, aes(x = lat_range_km, y = avg_distance_temp, color = scenario)) +
  geom_point(size = 0.1, alpha = 0.4) +
  geom_errorbar(aes(ymin = avg_distance_temp - se_distance_temp, 
                    ymax = avg_distance_temp + se_distance_temp),
                width = 0.1, alpha = 0.4) +
  scale_color_manual(values = c(
    "SSP126" = "#4DAF4A",
    "SSP245" = "#FFFF33",
    "SSP370" = "#FF7F00",
    "SSP585" = "#E41A1C" )) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(label=scales::comma,limits=c(0,10000),breaks=c(seq(0,10000,5000)),expand=c(0,0))+
  scale_y_continuous(label=scales::comma,limits = c(0,3200),breaks=(seq(0, 3000, 1000)),expand=c(0,0))+
  labs(x = "Range size (km)",
       y = "Distance to \ntemperature analog (km)") +
  theme_classic()

#all ssp as density plot

# Create the density plot with the combined dataframe
temp_analog <- ggplot(species_avg_all, aes(x = lat_range_km , y = avg_distance_temp, fill = scenario)) +
  stat_density_2d(geom = "polygon",bins = 50,contour = TRUE,alpha = 0.3) +
  scale_fill_manual(
    values = c(
      "SSP126" = "#4DAF4A",
      "SSP245" = "#FFFF33",
      "SSP370" = "#FF7F00",
      "SSP585" = "#E41A1C" ),
    name = "Scenario") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(label=scales::comma,limits=c(0,10000),breaks=c(seq(0,10000,5000)),expand=c(0,0))+
  scale_y_continuous(label=scales::comma,limits = c(0,2000),breaks=(seq(0, 2000, 1000)),expand=c(0,0))+
  labs(x = "Latitudinal range size (km)",
    y = "Distance to \ntemperature analog (km)") +
  theme_classic()


temp_analog <- ggplot(species_avg_all, aes(x = total_area/exceedall585$total_area[1]*100, y = avg_distance_temp, fill = scenario)) +
  stat_density_2d(geom = "polygon", bins = 50, contour = TRUE, alpha = 0.3) +
  scale_fill_manual(
    values = c(
      "SSP126" = "#4DAF4A",
      "SSP245" = "#FFFF33",
      "SSP370" = "#FF7F00",
      "SSP585" = "#E41A1C" ),
    name = "Scenario") +
  scale_x_continuous(label = scales::comma, limits = c(0, 70), breaks = c(seq(0, 70, 10)), expand = c(0, 0)) +
  scale_y_continuous(label = scales::comma, limits = c(0, 2000), breaks = (seq(0, 2000, 1000)), expand = c(0, 0)) +
  labs(x = "Range as proportion \nof domain area (%)",
       y = "Distance to \n temperature analog (km)") +
  theme_classic() +
  # theme(legend.position = c(0.8, 0.5),  # Position in the bottom right corner
  #       legend.background = element_rect(fill = "white", color = NA),
  #       legend.margin = margin(6, 6, 6, 6),
  #       legend.box.background = element_rect(color = "white"))+
  theme(legend.position = "none") +
  annotate("text", x = 60, y = 300, label = "SSP126", color = "black", fontface = "bold") +
  annotate("text", x = 60, y = 900, label = "SSP245", color = "black", fontface = "bold") +
  annotate("text", x = 60, y = 1300, label = "SSP370", color = "black", fontface = "bold") +
  annotate("text", x = 60, y = 1700, label = "SSP585", color = "black", fontface = "bold")







fig3_coral <- plot_grid(temp_lat,range_size, tempdelta_range,tempexceed_range,
                        pct_loss_compare_largest, temp_envelope,analog_loss,aragexceed_range,
                        labels = c( "A","B","C","D","E","F","G","H"), ncol = 4,nrow=2,align="hv")


ggsave(paste0(figures_folder, "fig3.pdf"), fig3_coral, 
       width = 15, height = 6, dpi = 300)  

ggsave(paste0(figures_folder, "fig3.png"), fig3_coral, 
       width = 15, height = 6, dpi = 350)

























