























scenarios <- c("ssp126","ssp245","ssp370","ssp585")


## for each ocean basin

#identify which species in which ocean
allcoral_coord <- fread( paste0(output_directory, 'allcoral_pixels_with_oceans_1982-1992.csv'))
species_coord$x_360 <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
ocean_reference <- unique(allcoral_coord[, .(x, y, ocean)])
setnames(ocean_reference, c("x", "y"), c("x_360", "y"))
species_coord_ocean <- merge(species_coord, ocean_reference, by=c("x_360", "y"), all.x=TRUE)

#find out where species are present
pacific_species <- unique(species_coord_ocean[ocean == "Pacific", .(species_id)])
atlantic_species <- unique(species_coord_ocean[ocean == "Atlantic", .(species_id)])
indian_species <- unique(species_coord_ocean[ocean == "Indian", .(species_id)])

indo_pacific <- merge(pacific_species, indian_species, by = "species_id")
pacific_atlantic <- merge(pacific_species, atlantic_species, by = "species_id") #empty
atlantic_indian <- merge(atlantic_species, indian_species, by = "species_id") #empty
all_oceans <- merge(pacific_atlantic, indian_species, by = "species_id") #empty

pacific_only <- pacific_species[!species_id %in% indo_pacific$species_id]
indian_only <- indian_species[!species_id %in% indo_pacific$species_id]


ocean_data <- fread(paste0(output_directory, "all_oceans_exceedance_summary_ssp585.csv"))
ocean_data[, pct_exceed_temp_area := temp_exceed_area / total_area]

# Extract allcoral data from ocean_data by ocean basin
allcoral_by_ocean <- ocean_data %>%
  filter(species_id == "allcoral") %>%
  select(ocean, year, scenario, pct_exceed_temp_area)

# Split into Atlantic, Pacific, and Indian
allcoral_atlantic <- allcoral_by_ocean %>% 
  filter(ocean == "Atlantic") %>%
  rename(allcoral_pct_temp = pct_exceed_temp_area)

allcoral_pacific <- allcoral_by_ocean %>% 
  filter(ocean == "Pacific") %>%
  rename(allcoral_pct_temp = pct_exceed_temp_area)

allcoral_indian <- allcoral_by_ocean %>% 
  filter(ocean == "Indian") %>%
  rename(allcoral_pct_temp = pct_exceed_temp_area)

allcoral_indopac <- allcoral_by_ocean %>% 
  filter(ocean == "indo-pacific") %>%
  rename(allcoral_pct_temp = pct_exceed_temp_area)

# Apply the comparison for each ocean basin separately
# For Atlantic
exceed585_atlantic <- exceed585 %>%
  # Filter to include only species found in Atlantic
  filter(species_id %in% atlantic_species$species_id) %>%
  # Join with Atlantic allcoral data
  left_join(
    allcoral_atlantic %>% select(year, scenario, allcoral_pct_temp),
    by = c("year", "scenario")
  ) %>%
  # Calculate the difference (species - allcoral)
  mutate(
    underestimated_pct_temp = pct_exceed_temp_area - allcoral_pct_temp,
    ocean = "Atlantic"
  ) %>%
  # Filter out any rows that might match allcoral
  filter(species_id != "allcoral")

# For Pacific
exceed585_pacific <- exceed585 %>%
  # Filter to include only species found in Pacific
  filter(species_id %in% pacific_only$species_id) %>%
  # Join with Pacific allcoral data
  left_join(
    allcoral_pacific %>% select(year, scenario, allcoral_pct_temp),
    by = c("year", "scenario")
  ) %>%
  # Calculate the difference (species - allcoral)
  mutate(
    underestimated_pct_temp = pct_exceed_temp_area - allcoral_pct_temp,
    ocean = "Pacific"
  ) %>%
  # Filter out any rows that might match allcoral
  filter(species_id != "allcoral")

# For Indian 
exceed585_indian <- exceed585 %>%
  # Filter to include only species found in Indian
  filter(species_id %in% indian_only$species_id) %>%
  # Join with Indian allcoral data
  left_join(
    allcoral_indian %>% select(year, scenario, allcoral_pct_temp),
    by = c("year", "scenario")
  ) %>%
  # Calculate the difference (species - allcoral)
  mutate(
    underestimated_pct_temp = pct_exceed_temp_area - allcoral_pct_temp,
    ocean = "Indian"
  ) %>%
  # Filter out any rows that might match allcoral
  filter(species_id != "allcoral")

# indo-pacific
exceed585_indopac <- exceed585 %>%
  # Filter to include only species found in indo-pacific
  filter(species_id %in% indo_pacific$species_id | 
           species_id %in% indian_only$species_id | 
           species_id %in% pacific_only$species_id) %>%
  # Join with indo-pacific allcoral data
  left_join(
    allcoral_indopac %>% select(year, scenario, allcoral_pct_temp),
    by = c("year", "scenario")
  ) %>%
  # Calculate the difference (species - allcoral)
  mutate(
    underestimated_pct_temp = pct_exceed_temp_area - allcoral_pct_temp,
    ocean = "indo-pacific"
  ) %>%
  # Filter out any rows that might match allcoral
  filter(species_id != "allcoral")

# Combine the datasets using rbind
exceed585_all_oceans <- rbind(exceed585_atlantic, exceed585_pacific, exceed585_indian,exceed585_indopac)

# Create a function to generate a plot for each ocean
create_ocean_plot <- function(data, ocean_name, show_x_label = FALSE) {
  # Filter data for the specific ocean
  ocean_data <- data %>% filter(ocean == ocean_name)
  
  # Calculate the mean underestimate for this ocean by year
  ocean_mean <- ocean_data %>%
    group_by(year) %>%
    summarize(mean_underestimate = mean(underestimated_pct_temp, na.rm = TRUE), .groups = "drop")
  
  # Create the plot
  ggplot() +
    # species difference lines
    geom_line(data = ocean_data,
              aes(x = year, 
                  y = underestimated_pct_temp*100,
                  group = species_id),
              color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
    # average of all species 
    geom_line(data = ocean_mean,
              aes(x = year, 
                  y = mean_underestimate*100),
              color = "blue",
              linewidth = 1.5) +
    # Zero reference line (where species = allcoral for this ocean)
    geom_hline(yintercept = 0, color = "red", linewidth = 1) +
    labs(x = ifelse(show_x_label, "Year", ""),
         y = "Temperature Area \nunderestimated (%)") +
    theme_classic() +
    theme(
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    scale_y_continuous(limits = c(-40, 100), 
                       breaks = seq(-40, 100, 20),
                       expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 5),
                       limits = c(2015, 2100))
}

# Create plots for each ocean
atlantic_plot <- create_ocean_plot(exceed585_all_oceans, "Atlantic")
pacific_plot <- create_ocean_plot(exceed585_all_oceans, "Pacific")
indian_plot <- create_ocean_plot(exceed585_all_oceans, "Indian", TRUE)
indopac_plot <- create_ocean_plot(exceed585_all_oceans, "indo-pacific", TRUE)




# habitat loss per ocean

create_habloss_ocean_plot <- function(ocean_name, exceed585_data, ocean_data) {
  # Filter data for the specific ocean
  ocean_species_data <- exceed585_data %>% 
    filter(species_id != "allcoral", ocean == ocean_name)
  
  ocean_summary_data <- ocean_species_data %>%
    group_by(year, scenario) %>%
    summarize(mean_unsuitable = mean(pct_exceed_temp_area*100, na.rm = TRUE), .groups = "drop")
  
  ocean_allcoral_data <- exceed585_data %>% 
    filter(species_id == "allcoral", ocean == ocean_name)
  
  ocean_line_data <- ocean_data %>%
    filter(ocean == ocean_name)
  
  # Create the plot
  ggplot() +
    # Species lines for this ocean
    geom_line(data = ocean_species_data,
              aes(x = year, y = pct_exceed_temp_area*100, group = interaction(species_id, scenario)),
              color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
    # Average of all species in this ocean
    geom_line(data = ocean_summary_data,
              aes(x = year, y = mean_unsuitable, group = scenario),
              color = "blue", linewidth = 1.5) +
    # Ocean-specific line
    geom_line(data = ocean_line_data,
              aes(x = year, y = pct_exceed_temp_area*100),
              color = "red", linewidth = 1.5) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0),
                       breaks = seq(0, 100, 20)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = "Year",
         y = "Habitat loss (area exceeding \ntemperature envelopes (%))",
         title = paste(stringr::str_to_title(ocean_name), "Ocean")) +
    theme_classic() +
    theme(
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", color = NA)
    )
}

# Create plots for each ocean
habloss_plot_atl <- create_habloss_ocean_plot("Atlantic", exceed585_all_oceans, ocean_data)
habloss_plot_pac <- create_habloss_ocean_plot("Pacific", exceed585_all_oceans, ocean_data)
habloss_plot_indopac <- create_habloss_ocean_plot("indo-pacific", exceed585_all_oceans, ocean_data)
habloss_plot_ind <- create_habloss_ocean_plot("Indian", exceed585_all_oceans, ocean_data)













#percentage of species above red line over time
species_above_threshold_temp_atl <- exceed585_atlantic %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_temp > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_temp_atl <- ggplot(species_above_threshold_temp_atl,
                                      aes(x = year, y = percent_above_threshold)) +
  geom_line(color = "blue", linewidth = 1.5) +
  labs(x = "Year",
       y = "Species \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))

species_above_threshold_temp_pac <- exceed585_pacific %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_temp > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_temp_pac <- ggplot(species_above_threshold_temp_pac,
                                      aes(x = year, y = percent_above_threshold)) +
  geom_line(color = "blue", linewidth = 1.5) +
  labs(x = "Year",
       y = "Species \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))


species_above_threshold_temp_ind <- exceed585_indian %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_temp > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_temp_ind <- ggplot(species_above_threshold_temp_ind,
                                      aes(x = year, y = percent_above_threshold)) +
  geom_line(color = "blue", linewidth = 1.5) +
  labs(x = "Year",
       y = "Species \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))


species_above_threshold_temp_indopac <- exceed585_indopac %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_temp > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_temp_indopac <- ggplot(species_above_threshold_temp_indopac,
                                      aes(x = year, y = percent_above_threshold)) +
  geom_line(color = "blue", linewidth = 1.5) +
  labs(x = "Year",
       y = "Species \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))


# percentage area of overall coral exceeded envelope vs percent area of species exceeded envelopes
temp_envelope_global <-  ggplot() +
  #species lines
  geom_line(data = exceed585 %>% filter(species_id != "allcoral"),
            aes(x = year, y = pct_exceed_temp_area*100, group = interaction(species_id, scenario)),
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_unsuitable = mean(pct_exceed_temp_area*100, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_unsuitable,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = pct_exceed_temp_area*100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100), expand=c(0,0),
                     breaks = seq(0, 100, 20)) +
  scale_x_continuous(expand=c(0,0))+
  labs(title="Global",x = "Year",
       y = "Habitat loss (area exceeding \ntemperature envelopes (%))") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )






fig_ocean <- plot_grid(habloss_plot_atl,atlantic_plot,pct_underestimated_temp_atl,
                       habloss_plot_indopac ,indopac_plot,pct_underestimated_temp_indopac,
                       temp_envelope_global,temp_envelope_underestimate,pct_underestimated_temp,
                        labels = c("A","B","C","D","E","F","G","H","I"),
                        ncol = 3, nrow=3, align="hv")

ggsave(paste0(figures_folder, "fig_S_temp_ocean.png"), fig_ocean, 
       width = 10, height = 10, dpi = 350)




# map of species that are overestimated 
library(rnaturalearth)
library(rasterVis)
library(sp)
library(gridExtra)


world_sf <- ne_coastline(scale = "medium", returnclass = "sf")

polygon <- st_polygon(x = list(rbind(c(-0.0001, 90),
                                     c(0, 90),
                                     c(0, -90),
                                     c(-0.0001, -90),
                                     c(-0.0001, 90)))) %>%
  st_sfc() %>%
  st_set_crs(4326)

world_sf_fix <- world_sf %>% st_difference(polygon)
world_sf_shifted <- st_shift_longitude(world_sf_fix)

# Convert sf to sp for use with sp.polygons
world_sp <- as(world_sf_shifted, "Spatial")


levelplot(overestimated_raster_atl, 
          main = "Atlantic",
          margin = list(FUN = mean), 
          col.regions = rev(topo.colors(100)),
          xlab = "Longitude", 
          ylab = "Latitude") + 
  latticeExtra::layer(sp::sp.polygons(world_sp, col = "black", lwd = 0.5, fill = NA))


overestimated_rast <- c(overestimated_raster_ind,overestimated_raster_pac,overestimated_raster_indopac, overestimated_raster_atl)
names(overestimated_rast) <- c("Indian", "Pacific", "IndoPacific", "Atlantic")

overestimate_all <- levelplot(overestimated_rast, 
          layout = c(2, 2), 
          col.regions = rev(topo.colors(100)),
          margin = list(
            FUN = mean,
            draw = TRUE,
            scales = list(x = list(draw = TRUE), y = list(draw = TRUE))
          ),
          xlab = "Longitude", 
          ylab = "Latitude") + 
  latticeExtra::layer(sp::sp.polygons(world_sp, col = "black", lwd = 0.5, fill = NA))




### species richness

species_richness <- species_coord_ocean %>%
  group_by(x_360, y) %>%
  summarize(richness = n_distinct(species_id), .groups = "drop")

species_richness$x_360 <- as.numeric(species_richness$x_360)
points_data <- vect(species_richness, geom = c("x_360", "y"), crs = "+proj=longlat +datum=WGS84")

richness_raster <- rasterize(points_data, r, field = "richness", fun = "max")

richness_map <- levelplot(richness_raster, 
          margin = list(FUN = mean), 
          col.regions = rev(topo.colors(100)),
          xlab = "Longitude", 
          ylab = "Latitude") + 
  latticeExtra::layer(sp::sp.polygons(world_sp, col = "black", lwd = 0.5, fill = NA))







#average range size 

species_with_range <- left_join(species_coord_ocean, exceed585_2100, by = "species_id")

avg_range_by_pixel <- species_with_range %>%
  group_by(x_360, y) %>%
  summarize(avg_range_size = mean(total_area, na.rm = TRUE),
            .groups = "drop")

avg_range_by_pixel$x_360 <- as.numeric(avg_range_by_pixel$x_360)
points_data <- vect(avg_range_by_pixel, geom = c("x_360", "y"), crs = "+proj=longlat +datum=WGS84")

range_raster <- rasterize(points_data, r, field = "avg_range_size", fun = "max")

range_map <- levelplot(range_raster, 
          margin = list(FUN = mean), 
          col.regions = rev(topo.colors(100)),
          xlab = "Longitude", 
          ylab = "Latitude") + 
  latticeExtra::layer(sp::sp.polygons(world_sp, col = "black", lwd = 0.5, fill = NA))







#scatterplot species richness x range size

richrange <- left_join(species_richness, avg_range_by_pixel, by = c("x_360", "y"))

richrange_plot <- ggplot(data = richrange, aes(x = avg_range_size, y = richness)) +
  #geom_point(alpha = 0.3) + 
  geom_pointdensity() + 
  # Use a blue gradient similar to your example
  scale_color_gradient(low = '#cbecff', high = "#2171b5", name = "Density") +
  labs(x = "Average range size\n (million km²)", 
       y = "Species richness") +
  scale_y_continuous(limits = c(-5, 600), expand=c(0,0),
                     breaks = seq(0, 600, 200)) +
  scale_x_continuous(expand=c(0,0),
                     limits=c(-5,25000000),
                     labels = function(x) format(x / 1000000, digits = 1, scientific = FALSE))+
  theme_classic() +
  theme(legend.position="none",
        aspect.ratio = 0.8)




species_summary <- species_with_range %>%
  group_by(species_id) %>%
  summarize(
    mean_latitude = mean(y, na.rm = TRUE),
    abs_latitude = abs(mean_latitude),  # Absolute latitude (distance from equator)
    total_range_size = sum(pixel_area, na.rm = TRUE),
    pct_habitat_loss = mean(pct_exceed_temp_area, na.rm = TRUE)
  )



p1 <- ggplot(species_summary, aes(x = abs_latitude, y = pct_habitat_loss*100)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "gray50", linetype = "dashed") +
  geom_point(aes(color = total_range_size/1000000), alpha = 0.6, size = 1.5) +
  # Custom color scale - small ranges in dark blue, large in light yellow
  scale_color_viridis(option = "plasma", direction = -1, name = "Range size\n (million km²)",labels = label_comma()) +
  labs(x = "Absolute latitude (degrees)",
    y = "Habitat loss (%)") +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12),
    legend.position = "none",
    aspect.ratio = 0.8)

# Residuals plot 
p2 <- ggplot(species_summary, aes(x = abs_latitude, y = residuals, color = total_range_size/1000000)) +
  # Add zero reference line
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  # Add points colored by range size
  geom_point(alpha = 0.6, size = 1.5) +
  # Custom color scale - small ranges in dark blue, large in light yellow
  scale_color_viridis(option = "plasma", direction = -1, name = "Range size",labels = label_comma()) +
  # Add conditional trend lines for small and large ranges
  geom_smooth(data = subset(species_summary, total_range_size < median(total_range_size)),
              aes(group = 1), method = "loess", color = "blue", size = 1.2) +
  labs(x = "Absolute latitude (degrees)",
    y = "Residual") +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 12),
    legend.position = "none" ,
    aspect.ratio = 0.8)

legend_plot <- ggplot(species_summary, aes(x = abs_latitude, y = pct_habitat_loss*100)) +
  geom_point(aes(color = total_range_size)) +
  scale_color_viridis(option = "plasma", direction = -1, 
                      name = "Range size (km²)", 
                      labels = label_comma()) +
  guides(color = guide_colorbar(barwidth = 10, barheight = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom")

# Extract just the legend
legend <- cowplot::get_legend(legend_plot)

range_lat <- ggplot(species_summary, aes(x = abs_latitude, y = total_range_size)) +
  geom_point( alpha = 0.3, size = 1.5) +
  # Use a log scale for range size
  scale_y_log10(labels = scales::comma) +
  labs(x = "Absolute latitude (degrees)",
    y = "Total range size (km²)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    aspect.ratio = 0.8)







#map of species with low loss (<85%)

species_with_low_loss <- species_with_range %>%
  filter(pct_exceed_temp_area < 0.85)

# Count species per pixel
species_count_by_pixel <- species_with_low_loss %>%
  group_by(x_360, y) %>%
  summarize(species_count = n(),
            .groups = "drop")

species_count_by_pixel$x_360 <- as.numeric(species_count_by_pixel$x_360)

points_data <- vect(species_count_by_pixel, geom = c("x_360", "y"), crs = "+proj=longlat +datum=WGS84")

count_raster <- rasterize(points_data, r, field = "species_count", fun = "max")
non_zero_mask <- count_raster > 0
count_raster_no_zeros <- mask(count_raster, non_zero_mask)



lowloss_map <- levelplot(count_raster,
                         main = "A",
                         margin = FALSE,#margin = list(FUN = mean), 
          col.regions = viridis(100, option = "heat"),  # Using the 'inferno' option from viridis
          xlab = "Longitude", 
          ylab = "Latitude") + 
  latticeExtra::layer(sp::sp.polygons(world_sp, col = "black", lwd = 0.5, fill = NA))



ggplot(species_with_low_loss, aes(x = y)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "black") +
  labs(x = "Latitude", y = "Count", title = "Distribution of Latitudes") +
  theme_classic()





# Convert richness raster to ggplot
richness_df <- as.data.frame(richness_raster, xy = TRUE)
colnames(richness_df)[3] <- "richness"

richness_ggplot <- ggplot() +
  geom_raster(data = richness_df, aes(x = x, y = y, fill = richness)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, label = scales::comma) +
  labs(fill = "Species \nrichness") +
  geom_sf(data = st_as_sf(world_sp), fill = NA, color = "grey", size = 0.5) +
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8))

# Convert range size raster to ggplot
range_df <- as.data.frame(range_raster, xy = TRUE)
colnames(range_df)[3] <- "avg_range_size"

range_ggplot <- ggplot() +
  geom_raster(data = range_df, aes(x = x, y = y, fill = avg_range_size/1000000)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, label = scales::comma) +
  labs(fill = "Average\nrange size\n(million km²)") +
  geom_sf(data = st_as_sf(world_sp), fill = NA, color = "grey", size = 0.5) +
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8))


# Convert the third map (low loss species count) to ggplot
lowloss_df <- as.data.frame(count_raster, xy = TRUE)
colnames(lowloss_df)[3] <- "species_count"

lowloss_ggplot <- ggplot() +
  geom_raster(data = lowloss_df, aes(x = x, y = y, fill = species_count)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, label = scales::comma) +
  labs(fill = "Species\ncount") +
  geom_sf(data = st_as_sf(world_sp), fill = NA, color = "grey", size = 0.5) +
  coord_sf() +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.8, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8))
    











library(ggpubr)

fig_s_maps <- ggarrange(
  # Top row
  ggarrange(richness_ggplot, range_ggplot, lowloss_ggplot, 
            labels = c("A", "B", "C"), 
            ncol = 3),
  
  # Bottom row
  ggarrange(
    # D and E
    ggarrange(richrange_plot, range_lat, 
              labels = c("D", "E"), 
              ncol = 2),
    
    # F and G with common legend
    ggarrange(p1, p2, 
              labels = c("F", "G"), 
              ncol = 2,
              common.legend = TRUE,
              legend = "bottom"),
    
    ncol = 2
  ),
  
  nrow = 2
)


# Save the combined plot
ggsave(paste0(figures_folder, "fig_s_maps.png"), fig_s_maps, 
       width = 15, height = 6, dpi = 350)


























