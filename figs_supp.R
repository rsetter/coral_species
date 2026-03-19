

### SUPP FIG Temperature


## C

exceed585_underestimate <- exceed585 %>%
  # Get the allcoral value for each year
  group_by(year, scenario) %>%
  mutate(
    # Get the allcoral percentage for comparison
    allcoral_pct_temp = pct_exceed_temp_area[species_id == "allcoral"],
    # Calculate the difference (species - allcoral)
    underestimated_pct_temp = pct_exceed_temp_area - allcoral_pct_temp,
    # pH
    allcoral_pct_ph = pct_exceed_ph_area[species_id == "allcoral"],
    underestimated_pct_ph = pct_exceed_ph_area - allcoral_pct_ph,
    #dual 
    allcoral_pct_dual = pct_exceed_dual_area[species_id == "allcoral"],
    underestimated_pct_dual = pct_exceed_dual_area - allcoral_pct_dual,
    allcoral_pct_dual_or = pct_dual_or_exceed_area[species_id == "allcoral"],
    underestimated_pct_dual_or = pct_dual_or_exceed_area - allcoral_pct_dual_or
  ) %>%
  ungroup() %>%
  # Filter out the allcoral row
  filter(species_id != "allcoral")

exceed585_underestimate %>%
  summarize(
    # Temperature
    temp_pct_species_affected = mean(underestimated_pct_temp > 0) * 100,
    temp_avg_magnitude = mean(underestimated_pct_temp[underestimated_pct_temp > 0], na.rm = TRUE) * 100,
    
    # pH
    ph_pct_species_affected = mean(underestimated_pct_ph > 0) * 100,
    ph_avg_magnitude = mean(underestimated_pct_ph[underestimated_pct_ph > 0], na.rm = TRUE) * 100,
    
    # Dual
    dual_pct_species_affected = mean(underestimated_pct_dual > 0) * 100,
    dual_avg_magnitude = mean(underestimated_pct_dual[underestimated_pct_dual > 0], na.rm = TRUE) * 100,
    dual_or_pct_species_affected = mean(underestimated_pct_dual_or > 0) * 100,
    dual_or_avg_magnitude = mean(underestimated_pct_dual_or[underestimated_pct_dual_or > 0], na.rm = TRUE) * 100
  )

temp_envelope_underestimate <- ggplot() +
  # species difference lines
  geom_line(data = exceed585_underestimate,
            aes(x = year, 
                y = underestimated_pct_temp*100,
                group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_underestimate %>%
              group_by(year) %>%
              summarize(mean_underestimate = mean(underestimated_pct_temp, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_underestimate*100),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold"))+
  scale_y_continuous(limits = c(-1, 100), 
                     breaks = seq(0, 100, 20),
                     expand=c(0,0)) +
  scale_x_continuous(expand=c(0,5),
                     limits=c(2015,2100))



#percentage of species above red line over time
species_above_threshold_temp <- exceed585_underestimate %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_temp > 0, na.rm = TRUE) * 100,
    .groups = "drop")




pct_underestimated_temp <- ggplot(species_above_threshold_temp,
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











# dhw asb year vs year
dhw585 <- rbind(dhwsp585,dhwall585)
dhw585$dhw_exceed_pct <- dhw585$dhw_exceed_area / dhw585$total_area
dhwsp585$dhw_exceed_pct <- dhwsp585$dhw_exceed_area / dhwsp585$total_area
dhwall585$dhw_exceed_pct <- dhwall585$dhw_exceed_area / dhwall585$total_area

dhw_asbarea <- ggplot() +
  geom_line(data=dhw585 %>% filter(species_id != "allcoral"),
            aes(x = year, y = dhw_exceed_pct*100, group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) + 
  # average of all species 
  geom_line(data = dhw585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_dhw_exceed = mean(dhw_exceed_pct*100, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_dhw_exceed),
            color = "blue",
            linewidth = 1.5) +
  geom_line(data=dhw585%>% filter(species_id == "allcoral"),
            aes(x = year, y = dhw_exceed_pct*100),
            color="red",linewidth=1.5)+
  theme_classic() +
  labs(x = "Year",
       y = "Area experiencing \nbleaching (%)") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))

dhw_asbdur <- ggplot() +
  geom_line(data=dhw585 %>% filter(species_id != "allcoral"),
            aes(x = year, y = mean_dhw_duration, group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) + 
  # average of all species 
  geom_line(data = dhw585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_dhw_duration = mean(mean_dhw_duration, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_dhw_duration),
            color = "blue",
            linewidth = 1.5) +
  geom_line(data=dhw585%>% filter(species_id == "allcoral"),
            aes(x = year, y = mean_dhw_duration),
            color="red",linewidth=1.5)+
  theme_classic() +
  labs(x = "Year",
       y = "Mean duration \n(months)") +
  theme( legend.position = "right",
         panel.grid.minor = element_blank())+
  scale_y_continuous(limits = c(-0.1, 12), 
                     breaks = seq(0, 12, 2),
                     expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 5),
                     limits = c(2015, 2100))








fig_s_temp_dhw <- plot_grid(temp_envelope_underestimate, pct_underestimated_temp,
                            atlantic_plot,pct_underestimated_temp_atl,
                            pacific_plot,pct_underestimated_temp_pac,
                            indian_plot,pct_underestimated_temp_ind,
                            dhw_asbarea,dhw_asbdur,
                            labels = c( "A","B","C","D","E","F","G","H","I","J"),
                            ncol = 2,nrow=5, align="hv")


ggsave(paste0(figures_folder, "fig_S_temp_dhw.png"), fig_s_temp_dhw, 
       width = 10, height = 15, dpi = 350)












#### SUPP FIG pH



### D

# percentage area of overall coral exceeded envelope vs percent area of species exceeded envelopes


ph_envelope <-  ggplot() +
  #species lines
  geom_line(data = exceed585 %>% filter(species_id != "allcoral"),
            aes(x = year, 
                y = pct_exceed_ph_area*100,
                group = interaction(species_id, scenario)),
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(pct_unsuitable = mean(pct_exceed_ph_area, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = pct_unsuitable*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = pct_exceed_ph_area*100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_viridis_d() +
  labs(x = "Year",
       y = "Area below \npH envelopes (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )



### E

exceedall585$mean_ph <- exceedall585$mean_ph_diff + exceedall585$ph_min
exceedall585$mean_ph_delta <- exceedall585$mean_ph - exceedall585$ph_mean
exceedsp585$mean_ph <- exceedsp585$mean_ph_diff + exceedsp585$ph_min
exceedsp585$mean_ph_delta <- exceedsp585$mean_ph - exceedsp585$ph_summer_mean

ph_delta <-  ggplot() +
  #species lines
  geom_line(data = exceedsp585 %>% filter(species_id != "allcoral"),
            aes(x = year, 
                y = mean_ph_delta,
                group = interaction(species_id, scenario)),
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceedsp585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_summer_delta = mean(mean_ph_delta, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_summer_delta),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceedall585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = mean_ph_delta,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(-0.45, 0.0)) +
  labs(x = "Year",
       y = "Summer pH \ndifference") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )



## F

ph_envelope_underestimate <- ggplot() +
  # species difference lines
  geom_line(data = exceed585_underestimate,
            aes(x = year, 
                y = underestimated_pct_ph*100,
                group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_underestimate %>%
              group_by(year) %>%
              summarize(mean_underestimate = mean(underestimated_pct_ph, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_underestimate*100),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))+
  scale_y_continuous(limits = c(-0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) 

#percentage of species above red line over time
species_above_threshold_ph <- exceed585_underestimate %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_ph > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_ph <- ggplot(species_above_threshold_ph,
                                aes(x = year, y = percent_above_threshold)) +
  geom_line(color = "darkblue", linewidth = 1.5) +
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

## ph profiles
ph_lat <- ggplot(combined_lat585, aes(x = latitude, y = ph, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Mean summer \npH",
       color = "Scenario") +
  scale_x_continuous(limits=c(-40,45))
#theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#ph increase v range size
phdelta_range <- ggplot(exceedsp585[exceedsp585$year==2100,], aes(x = total_area/exceedall585$total_area[1]*100, y = mean_ph_delta)) +
  geom_pointdensity()+
  scale_color_gradient(low='#cbecff',high="#2171b5") + 
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain size (%)", 
       y = "Mean ph \ndecrease")+
  scale_x_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20))+
  theme(legend.position="none")


# % loss vs range size
phexceed_range <- ggplot(exceed585[exceed585$year==2100,], aes(x = total_area/exceedall585$total_area[1]*100, y = pct_exceed_ph_area*100)) +
  geom_pointdensity()+
  scale_color_gradient(low='#cbecff',high="#2171b5") + #scale_color for point, scale_fill for hex
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Habitat loss \n(% area below today's min pH)")+
  scale_x_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20))+
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20))+
  theme(legend.position="none")




ph_analog <- ggplot(species_avg, aes(x = total_area, y = avg_distance_ph)) +
  geom_point(size = 0.1, alpha = 0.2) +  # Small points with high transparency
  geom_errorbar(aes(ymin = avg_distance_ph - se_distance_ph, 
                    ymax = avg_distance_ph + se_distance_ph),
                width = 0.1, alpha = 0.2) +  # Error bars with transparency
  labs(x = "Range size (km)",
       y = "Distance to \npH analog (km)") +
  theme_classic() +
  scale_x_continuous(label=scales::comma,limits=c(0,20000),breaks=c(seq(0,20000,5000)),expand=c(0,500))+
  scale_y_continuous(label=scales::comma,limits = c(0,25000),breaks=(seq(0, 25000, 5000)),expand=c(0,0))



fig_s_ph <- plot_grid( ph_lat,ph_delta, phdelta_range,phexceed_range,
                       ph_envelope, ph_envelope_underestimate, pct_underestimated_ph, 
                       labels = c( "A","B","C","D","E","F","G"),
                       ncol = 2,nrow=4, align="hv")


ggsave(paste0(figures_folder, "fig_S_ph.png"), fig_s_ph, 
       width = 8, height = 9, dpi = 350)










### SUPP FIG arag
exceedsp126_arag <- fread(paste0(output_directory, "species_exceedance_summary_arag_", scenarios[[1]], ".csv"))
exceedsp245_arag <- fread(paste0(output_directory, "species_exceedance_summary_arag_", scenarios[[2]], ".csv"))
exceedsp370_arag <- fread(paste0(output_directory, "species_exceedance_summary_arag_", scenarios[[3]], ".csv"))
exceedsp585_arag <- fread(paste0(output_directory, "species_exceedance_summary_arag_", scenarios[[4]], ".csv"))
exceedall126_arag <- fread(paste0(output_directory, "allcoral_exceedance_summary_arag_", scenarios[[1]], ".csv"))
exceedall245_arag <- fread(paste0(output_directory, "allcoral_exceedance_summary_arag_", scenarios[[2]], ".csv"))
exceedall370_arag <- fread(paste0(output_directory, "allcoral_exceedance_summary_arag_", scenarios[[3]], ".csv"))
exceedall585_arag <- fread(paste0(output_directory, "allcoral_exceedance_summary_arag_", scenarios[[4]], ".csv"))

species_envelope_arag <- fread(paste0(output_directory, 'climate_envelopes_1982-1992_with_arag.csv'))
species_envelope_arag$id_no <- as.character(species_envelope_arag$id_no)
exceedsp585_arag$species_id <- as.character(exceedsp585_arag$species_id)
exceedsp585_arag <- merge(exceedsp585_arag,species_envelope_arag,by.x="species_id",by.y="id_no",all.x=T)
allcoral_envelope_arag <- fread(paste0(output_directory, 'allcoral_envelopes_1982-1992_with_arag.csv'))
exceedall585_arag <- merge(exceedall585_arag,allcoral_envelope_arag,by.x="species_id",by.y="id_no",all.x=T)
exceed585_arag <- rbind(exceedsp585_arag,exceedall585_arag,fill=T)
exceed585_arag$pct_exceed_arag_area <- exceed585_arag$arag_exceed_area  / exceed585_arag$total_area


arag_envelope <-  ggplot() +
  #species lines
  geom_line(data = exceed585_arag %>% filter(species_id != "allcoral"),
            aes(x = decade, 
                y = pct_exceed_arag_area*100,
                group = interaction(species_id, scenario)),
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_arag %>% 
              filter(species_id != "allcoral") %>%
              group_by(decade, scenario) %>%
              summarize(pct_unsuitable = mean(pct_exceed_arag_area, na.rm = TRUE), .groups = "drop"),
            aes(x = decade, 
                y = pct_unsuitable*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585_arag %>% filter(species_id == "allcoral"),
            aes(x = decade, 
                y = pct_exceed_arag_area*100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_viridis_d() +
  labs(x = "Year",
       y = "Area below \nΩarag envelopes (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )




exceed585_arag$mean_arag <- exceed585_arag$mean_arag_diff + exceed585_arag$aragonite_min 
exceed585_arag$mean_arag_delta <- exceed585_arag$mean_arag - exceed585_arag$aragonite_mean 

arag_delta <-  ggplot() +
  #species lines
  geom_line(data = exceed585_arag %>% filter(species_id != "allcoral"),
            aes(x = decade, 
                y = mean_arag_delta,
                group = interaction(species_id, scenario)),
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_arag %>% 
              filter(species_id != "allcoral") %>%
              group_by(decade, scenario) %>%
              summarize(mean_summer_delta = mean(mean_arag_delta, na.rm = TRUE), .groups = "drop"),
            aes(x = decade, 
                y = mean_summer_delta),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585_arag %>% filter(species_id == "allcoral"),
            aes(x = decade, 
                y = mean_arag_delta,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  #scale_y_continuous(limits = c(-0.45, 0.0)) +
  labs(x = "Year",
       y = "Ωarag difference") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )



exceed585_underestimate_arag <- exceed585_arag %>%
  # Get the allcoral value for each year
  group_by(decade, scenario) %>%
  mutate(
    # Get the allcoral percentage for comparison
    allcoral_pct_arag = pct_exceed_arag_area [species_id == "allcoral"],
    # Calculate the difference (species - allcoral)
    underestimated_pct_arag = pct_exceed_arag_area - allcoral_pct_arag,
  ) %>%
  ungroup() %>%
  # Filter out the allcoral row
  filter(species_id != "allcoral")

exceed585_underestimate_arag %>%
  summarize(
    # Temperature
    temp_pct_species_affected = mean(underestimated_pct_arag > 0) * 100,
    temp_avg_magnitude = mean(underestimated_pct_arag[underestimated_pct_arag > 0], na.rm = TRUE) * 100,
  )


arag_envelope_underestimate <- ggplot() +
  # species difference lines
  geom_line(data = exceed585_underestimate_arag,
            aes(x = decade, 
                y = underestimated_pct_arag*100,
                group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_underestimate_arag %>%
              group_by(decade) %>%
              summarize(mean_underestimate = mean(underestimated_pct_arag, na.rm = TRUE), .groups = "drop"),
            aes(x = decade, 
                y = mean_underestimate*100),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area \nunderestimated (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))+
  scale_y_continuous(limits = c(-0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) 

#percentage of species above red line over time
species_above_threshold_arag <- exceed585_underestimate_arag %>%
  group_by(decade) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_arag > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_arag <- ggplot(species_above_threshold_arag,
                                  aes(x = decade, y = percent_above_threshold)) +
  geom_line(color = "darkblue", linewidth = 1.5) +
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




calculate_aragonite_profile <- function(output_directory, scenario, cmip_folder) {
  # Read current coral data
  allcoral_coord <- fread(paste0(output_directory, 'allcoral_pixels_1982-1992.csv'))
  
  # Read historical aragonite data and extract layer 16 (current)
  arag_hist <- rast(paste0(cmip_folder, "arag/Aragonite_median_historical_foc.nc"))
  arag_current <- arag_hist[[16]]
  
  # Read future aragonite data and extract the last layer (2100)
  arag_future <- rast(paste0(cmip_folder, "arag/Aragonite_median_", scenario, "_foc.nc"))
  arag_2100 <- arag_future[[9]]  # Last layer (9th decade)
  
  # Determine hemisphere for each coordinate
  allcoral_coord[, hemisphere := ifelse(y > 0, "N", "S")]
  
  # Extract current aragonite values
  arag_current_values <- terra::extract(arag_current, allcoral_coord[, .(x, y)])
  allcoral_coord[, arag_current := arag_current_values[,2]]
  
  # Extract future aragonite values
  arag_future_values <- terra::extract(arag_2100, allcoral_coord[, .(x, y)])
  allcoral_coord[, arag_future := arag_future_values[,2]]
  
  # Calculate current latitude data
  current_lat <- allcoral_coord[, .(
    aragonite = mean(arag_current, na.rm = TRUE)
  ), by = y]
  setnames(current_lat, "y", "latitude")
  current_lat[, scenario := "Today"]
  
  # Calculate future latitude data
  future_lat <- allcoral_coord[, .(
    aragonite = mean(arag_future, na.rm = TRUE)
  ), by = y]
  setnames(future_lat, "y", "latitude")
  future_lat[, scenario := "Future"]
  
  # Combine current and future data
  combined_lat <- rbind(current_lat, future_lat)
  
  return(combined_lat)
}

# Apply the function for each scenario
for(scenario in scenarios){
  var_name <- paste0('combined_lat', gsub('ssp','', scenario), '_arag')
  assign(var_name, calculate_aragonite_profile(output_directory, scenario, cmip_folder), envir = .GlobalEnv)
}





## aragonite profiles
arag_lat <- ggplot(combined_lat585_arag, aes(x = latitude, y = aragonite, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Mean Ωarag",
       color = "Scenario") +
  scale_x_continuous(limits=c(-40,45))
#theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())



#aragonite increase v range size
aragdelta_range <- ggplot(exceed585_arag[exceed585_arag$decade==2100,], aes(x = total_area/exceedall585_arag$total_area[1]*100, y = mean_arag_delta)) +
  geom_pointdensity()+
  scale_color_gradient(low='#cbecff',high="#2171b5") + 
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain size (%)", 
       y = "Mean aragonite \ndecrease")+
  scale_x_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20))+
  theme(legend.position="none")




# Sort the Aragonite data by range size (as a proportion of domain)
sorted_dataarag <- exceed585_arag[exceed585_arag$decade == 2100, ] %>%
  filter(species_id != "allcoral") %>% 
  arrange(total_area / exceedall585_arag$total_area[1])

# Calculate cumulative percentage of species
sorted_dataarag$cum_species_pct <- seq_len(nrow(sorted_dataarag)) / nrow(sorted_dataarag)

# % loss vs range size
aragexceed_range <- ggplot(exceed585_arag[exceed585_arag$decade==2100,] %>% filter(species_id != "allcoral"), 
                           aes(x = total_area/exceedall585_arag$total_area[1]*100, y = pct_exceed_arag_area*100)) +
  geom_pointdensity()+
  scale_color_gradient(low='#cbecff',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  # Add the cumulative species line
  geom_line(data = sorted_dataarag, aes(x = total_area/exceedall585_arag$total_area[1]*100, y = cum_species_pct*100), color = "blue", linewidth = 1) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain size (%)", 
       y = "Area exceeding \nΩarag envelope (%)")+
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),expand = c(0,0))+
  scale_y_continuous(
    name = "Area exceeding\nΩarag envelope (%)",
    sec.axis = sec_axis(~., name = "Cumulative species (%)"),
    expand=c(0,0),
    limits = c(0, 100.5),
    breaks = seq(0, 100, 20)
  ) +
  theme(legend.position="none")+
  theme(legend.position="none",
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.line.y.right = element_line(color="blue"))






#this hasn't been calculated yet
# arag_analog <- ggplot(species_avg_arag, aes(x = range_km, y = avg_distance_arag)) +
#   geom_point(size = 0.1, alpha = 0.2) +  # Small points with high transparency
#   geom_errorbar(aes(ymin = avg_distance_arag - se_distance_arag, 
#                     ymax = avg_distance_arag + se_distance_arag),
#                 width = 0.1, alpha = 0.2) +  # Error bars with transparency
#   labs(x = "Range size (km)",
#        y = "Distance to \naragonite analog (km)") +
#   theme_classic() +
#   scale_x_continuous(label=scales::comma,limits=c(0,20000),breaks=c(seq(0,20000,5000)),expand=c(0,500))+
#   scale_y_continuous(label=scales::comma,limits = c(0,25000),breaks=(seq(0, 25000, 5000)),expand=c(0,0))



fig_s_arag <- plot_grid(arag_lat, arag_delta,aragdelta_range, aragexceed_range,
                        arag_envelope, arag_envelope_underestimate, pct_underestimated_arag, 
                        labels = c("A","B","C","D","E","F","G"),
                        ncol = 2, nrow=4, align="hv")


ggsave(paste0(figures_folder, "fig_S_arag.png"), fig_s_arag, 
       width = 8, height = 9, dpi = 350)








#dataframe with species that have <80% habitat loss for either temp or arag
low_impact_species <- data.frame()
arag_data <- exceed585_arag[exceed585_arag$decade==2100,]
temp_data <- exceed585[exceed585$year==2100,]
temparag_data <- merge(arag_data, temp_data, by="species_id", suffixes=c("_arag", "_temp"))

low_impact_species <- temparag_data[temparag_data$pct_exceed_arag_area*100 < 80 | 
                                      temparag_data$pct_exceed_temp_area*100 < 80, ]

ggplot(low_impact_species, aes(x = pct_exceed_temp_area*100, y = pct_exceed_arag_area*100)) +
  geom_point() +
  labs(x = "Temperature habitat loss (%)",
       y = "Aragonite habitat loss (%)") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  theme_classic() +
  theme(legend.position = "none")


ggplot(temparag_data, aes(x = pct_exceed_temp_area*100, y = pct_exceed_arag_area*100)) +
  geom_pointdensity() +
  scale_color_gradient(low = '#cbecff', high = "#2171b5") +
  # Add labels and formatting
  labs(x = "Temperature habitat loss (%)",
       y = "Aragonite habitat loss (%)") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  theme_classic() +
  theme(legend.position = "none")




### SUPP FIG dual



## G

#either temperature or ph envelopes exceeded 
dual_envelope <-  ggplot() +
  #species lines
  geom_line(data = exceed585 %>% filter(species_id != "allcoral"),
            aes(x = year, 
                y = pct_dual_or_exceed_area*100,
                group = interaction(species_id, scenario)), 
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_envelope = mean(pct_dual_or_exceed_area, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_envelope*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = pct_dual_or_exceed_area*100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20)) +  scale_color_viridis_d() +
  labs(x = "Year",
       y = "Area exceeding \ntemperature and pH envelopes (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))





#H 

dual_envelope_underestimate <- ggplot() +
  # species difference lines
  geom_line(data = exceed585_underestimate,
            aes(x = year, 
                y = underestimated_pct_dual*100,
                group = species_id),
            color = "#4292c6",alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585_underestimate %>%
              group_by(year) %>%
              summarize(mean_underestimate = mean(underestimated_pct_dual, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_underestimate*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area underestimated by \nhabitat-level dual analysis (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))+
  scale_y_continuous(limits = c(-0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) 




#percentage of species above red line over time
species_above_threshold_dual <- exceed585_underestimate %>%
  group_by(year) %>%
  summarize(
    percent_above_threshold = mean(underestimated_pct_dual_or  > 0, na.rm = TRUE) * 100,
    .groups = "drop")

pct_underestimated_dual <- ggplot(species_above_threshold_dual,
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


fig_s_dual <- plot_grid( dual_envelope, dual_envelope_underestimate, pct_underestimated_dual,
                         labels = c( "A","B","C"),
                         ncol = 3,nrow=1, align="hv")


ggsave(paste0(figures_folder, "fig_S_dual.png"), fig_s_dual, 
       width = 9, height = 3, dpi = 350)














combined_lathist <- combined_lat585 %>% 
  filter(scenario == "Today") %>% 
  mutate(scenario = "Historical",
         ssp="Historical")

# Combine with future scenarios from each dataset
combined_latssp <- rbind(
  combined_lathist,
  combined_lat126 %>% filter(scenario == "Future") %>% mutate(ssp = "SSP126"),
  combined_lat245 %>% filter(scenario == "Future") %>% mutate(ssp = "SSP245"),
  combined_lat370 %>% filter(scenario == "Future") %>% mutate(ssp = "SSP370"),
  combined_lat585 %>% filter(scenario == "Future") %>% mutate(ssp = "SSP585"))

temp_lat_ssp <- ggplot(combined_latssp, aes(x = latitude, y = temperature, color = ssp)) +
  geom_line(size = 1.2) +
  scale_color_manual(
    values = c(
      "Historical" = "black",
      "SSP126" = "#4DAF4A",  
      "SSP245" = "#f9ec4e",  
      "SSP370" = "#FF7F00", 
      "SSP585" = "#E41A1C"  
    ),
    name = "Scenario"
  ) +
  theme_classic() +
  labs(x = "Latitude (degrees)", 
       y = "Mean summer \ntemperature (°C)") +
  scale_x_continuous(limits = c(-45, 45), expand = c(0, 0), breaks = seq(-45, 45, 45)) +
  scale_y_continuous(limits = c(15, 35), breaks = seq(15, 35, 5), expand = c(0, 0)) +
  theme(
    legend.position = c(0.5, 0.1),
    legend.justification = c(0.5, 0),
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5),
    legend.margin = margin(6, 6, 6, 6)
  )






#critical range threshold 

displacement585 <- calculate_latitudinal_displacement(combined_lat585 %>% filter(scenario == "Today") %>% select(latitude,temperature),
                                                   combined_lat585 %>% filter(scenario == "Future") %>% select(latitude,temperature))
# Filter out cross-hemisphere matches
displacement585 <- displacement585 %>%
  mutate(cross_hemisphere = !is.na(latitude) & !is.na(future_lat) & sign(latitude) != sign(future_lat)) %>%
  # Set displacement to NA for cross-hemisphere matches
  mutate(displacement = ifelse(cross_hemisphere, NA, displacement),
         display_displacement = ifelse(cross_hemisphere, NA, display_displacement))

min_critical_range <- min(displacement585$display_displacement, na.rm = TRUE)
max_critical_range <- max(displacement585$display_displacement, na.rm = TRUE)


# Plot displacement vs. latitude
ggplot(displacement585, aes(x = latitude, y = displacement)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Displacement (degrees latitude)",
       title = "Required latitudinal shift to maintain current temperature")+
  scale_x_continuous(limits = c(-50, 50))+
  scale_y_continuous(limits=c(0,30))



# latitude and distance to analog (critical range threshold)



scenarios <- c("ssp126", "ssp245", "ssp370", "ssp585")

analog_latitude <- function(scenario, output_directory_analog) {
  # Load North and South data
  diffanalogN <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerN', "_", scenario, "_rangeTemp", 1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
  diffanalogS <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeanfocsummerS', "_", scenario, "_rangeTemp", 1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
  
  # Take first observation for each coordinate pair
  diffanalogN <- diffanalogN[, .SD[1], by = .(px_temp, py_temp)]
  diffanalogS <- diffanalogS[, .SD[1], by = .(px_temp, py_temp)]
  
  # Merge with coral data - using your approach with x,y instead of px_temp,py_temp
  north_result <- merge(north_coral, diffanalogN, by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
  south_result <- merge(south_coral, diffanalogS, by.x = c("x", "y"), by.y = c("px_temp", "py_temp"), all.x = TRUE)
  
  # Combine results and add metadata
  all_result <- rbind(north_result, south_result)
  all_result[, scenario := scenario]
  all_result[, latitude := round(y)]
  
  # Calculate summary statistics by latitude
  lat_summary <- all_result[!is.na(distance_temp), .(
    avg_distance = mean(distance_temp, na.rm = TRUE),
    median_distance = median(distance_temp, na.rm = TRUE),
    sd_distance = sd(distance_temp, na.rm = TRUE),
    min_distance = min(distance_temp, na.rm = TRUE),
    max_distance = max(distance_temp, na.rm = TRUE),
    n_points = .N
  ), by = .(latitude, scenario)]
  
  return(list(full_data = all_result, summary = lat_summary))
}


north_coral <- allcoral_coord[y >= 0]
south_coral <- allcoral_coord[y < 0]

ssp_analoglat <- list()
ssp_analoglat_summary <- list()
for (scenario in scenarios) {
  results <- analog_latitude(scenario, output_directory_analog)
  ssp_analoglat[[scenario]] <- results$full_data
  ssp_analoglat_summary[[scenario]] <- results$summary
}
ssp_analoglat_summary <- rbindlist(ssp_analoglat_summary)

sspanaloglat <- ggplot(ssp_analoglat_summary, aes(x = latitude, y = avg_distance, color = scenario)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c("ssp126" = "#4DAF4A", "ssp245" = "#f9ec4e", "ssp370" = "#FF7F00", "ssp585" = "#E41A1C"), name = "Scenario") +
  labs(x = "Latitude (degrees)", y = "Distance to analog (km)") +
  scale_x_continuous(limits = c(-50, 50))+
  scale_y_continuous(limits=c(0,2000))+
  theme_classic()+
  theme(legend.position = "none")


fig_s_loss_ssp <- plot_grid( temp_lat_ssp, sspanaloglat, 
                             labels = c( "A","B"),
                             ncol = 2,nrow=1, align="hv")

ggsave(paste0(figures_folder, "fig_S_loss_ssp.png"), fig_s_loss_ssp, 
       width = 9, height = 4, dpi = 350)





















#habitat loss by ssp
exceedsp126$scenario <- "SSP126"
exceedsp245$scenario <- "SSP245"
exceedsp370$scenario <- "SSP370"
exceedsp585$scenario <- "SSP585"

exceed_ssp <- rbind(
  exceedsp126[exceedsp126$year == 2100, ],
  exceedsp245[exceedsp245$year == 2100, ],
  exceedsp370[exceedsp370$year == 2100, ],
  exceedsp585[exceedsp585$year == 2100, ])
exceed_ssp$pct_exceed_temp_area <- exceed_ssp$temp_exceed_area / exceed_ssp$total_area

ref_area <- exceedall585$total_area[1]

loss_ssp <- ggplot(exceed_ssp, aes(x = total_area/ref_area*100, 
                       y = pct_exceed_temp_area*100, 
                       fill = scenario)) +
  # Use stat_density_2d with polygon geometry for density contours
  stat_density_2d(geom = "polygon", 
                  bins = 9, 
                  contour = TRUE, 
                  alpha = 0.3, 
                  h = c(12,12),
                  position = "identity" ,
                  aes(fill = scenario)) +
  scale_fill_manual(
    values = c(
      "SSP126" = "#4DAF4A",
      "SSP245" = "#f9ec4e",
      "SSP370" = "#FF7F00",
      "SSP585" = "#E41A1C" ),
    name = "Scenario") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand = c(0, 0)) +
  labs(x = "Range as proportion\n of domain size (%)",
       y = "Area exceeding \nenvelope (%)") +
  theme_classic() +
  theme(legend.position = "right")











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
  scale_fill_virifdis_c(name = "Number of\nspecies",
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
