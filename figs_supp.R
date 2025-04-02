







### FIG S_?

#percentage of area per species that remains suitable refugia over time and environmental conditions at those sites
#Quantification of how different future climate will be - a measure of “adaptation potential” or “resilience requirement”

#process species files
adapt585 <- list.files(paste0(output_directory,"threshold/"), 
                       pattern = "species_.*_exceedance_ssp585.csv",
                       full.names = TRUE) %>%
  map_df(read_csv, show_col_types = FALSE) %>%
  group_by(species_id, year) %>%
  summarise(
    
    total_area = sum(pixel_area,na.rm=T),
    
    #temp
    area_suitable_temp = sum(pixel_area[temp_exceed == 0],na.rm=T),
    pct_suitable_temp = (area_suitable_temp / total_area) * 100,
    suittemp_exceedph = mean(ph_exceed[temp_exceed == 0], na.rm = T),
    
    #ph
    area_suitable_ph = sum(pixel_area[ph_exceed == 0],na.rm=T),
    pct_suitable_ph = (area_suitable_ph / total_area) * 100,
    suitph_exceedtemp = mean(temp_exceed[ph_exceed == 0], na.rm = T),
    .groups = "drop"
  )



adapt_ph <- ggplot(adapt585, aes(x = year, y = pct_suitable_temp, group = species_id)) +
  geom_line(aes(color = suittemp_exceedph),alpha = 0.6) +
  scale_color_gradientn(colors = c("red","orange","yellow","grey80"),
                        name = "Mean pH\nexceedance") +
  labs(x = "Year", y = "% Area within temperature envelope") +
  theme_bw()

adapt_temp <- ggplot(adapt585, aes(x = year, y = pct_suitable_ph, group = species_id)) +
  geom_line(aes(color = suitph_exceedtemp),alpha = 0.6) +
  scale_color_gradientn(colors = c("grey80","yellow","orange","red"),
                        name = "Mean temperature\nexceedance") +
  #  scale_fill_distiller(palette = "YlGnBu", direction = 1, label=scales::comma) +
  labs(x = "Year", y = "% Area within pH envelope") +
  theme_bw()

exceed_temp <- ggplot(exceed585, aes(x = year, y = pct_exceed_temp_area, group = species_id)) +
  geom_line(aes(color = mean_temp_exceed),alpha = 0.6) +
  scale_color_gradientn(colors = c("grey80","yellow","orange","red"),
                        name = "Mean temperature\nexceedance") +
  labs(x = "Year", y = "% area exceed temperature envelope") +
  theme_bw()

exceed_ph <- ggplot(exceed585, aes(x = year, y = pct_exceed_ph_area, group = species_id)) +
  geom_line(aes(color = mean_ph_exceed),alpha = 0.6) +
  scale_color_gradientn(colors = c("red","orange","yellow","grey80"),
                        name = "Mean pH\nexceedance") +
  labs(x = "Year", y = "% area exceed pH envelope") +
  theme_bw()


plot_grid(exceed_temp, exceed_ph, labels = c("A", "B"), ncol = 2)






### FIG 3

# map - ssp585 2100 number of species exceeding envelopes
# hotspots of of multiple species exceeding thresholds

#open all species pixel-level ssp585 2100 
species_files <- list.files(path = file.path(output_directory, "threshold"), 
                            pattern = paste0("_exceedance_", "ssp585", ".csv"),
                            full.names = TRUE)
species_files <- species_files[-c(1)]
all_species_data <- rbindlist(lapply(species_files, function(file) {
  fread(file)[year == 2100]
}))
#remove NA rows
all_species_data_valid <- all_species_data[!is.na(all_species_data$temp_exceed),]

# aggregate counts by pixel
pixel_counts <- all_species_data_valid[, .(
  temp_count = sum(temp_exceed > 0, na.rm = TRUE),
  ph_count = sum(ph_exceed < 0, na.rm = TRUE),
  dual_count = sum(dual_exceed == TRUE, na.rm = TRUE),
  total_species = .N #how many species in a given pixel
), by = .(x, y)]

pixel_counts[, `:=`(
  temp_percent = 100 * temp_count / total_species,
  ph_percent = 100 * ph_count / total_species,
  dual_percent = 100 * dual_count / total_species
)]

ggplot()+
  geom_point(data=pixel_counts,aes(x=total_species,y=dual_count))


world <- map_data('world',wrap=c(0,360))

num_dual_map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "darkgrey", color = "darkgrey", size = 0.2) +
  geom_tile(data = pixel_counts, 
            aes(x = x, y = y, fill = dual_count)) +
  scale_fill_viridis_c(name = "Number of\nspecies",
                       option = "magma",
                       direction = -1) +
  coord_fixed(expand = FALSE) +
  theme_void()+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey85", color = NA),
        plot.title = element_text(hjust = 0.5))

pct_dual_map <- ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "darkgrey", color = "darkgrey", size = 0.2) +
  geom_tile(data = pixel_counts, 
            aes(x = x, y = y, fill = dual_percent)) +
  scale_fill_viridis_c(name = "Percent of\nspecies",
                       option = "magma",
                       direction = -1) +
  coord_fixed(expand = FALSE) +
  theme_void() +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey85", color = NA),
        plot.title = element_text(hjust = 0.5))

plot_grid(num_dual_map, pct_dual_map, labels = c("A", "B"), ncol = 2)

species_dual_exceed <- all_species_data %>%
  filter(dual_exceed == TRUE) %>%
  distinct(species_id) %>%
  nrow()
