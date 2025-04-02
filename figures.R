# figures
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggridges)


source("functions.R")

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
exceed585$pct_exceed1_ph_area <- exceed585$ph_exceed1_area / exceed585$total_area
exceed585$pct_exceed2_ph_area <- exceed585$ph_exceed2_area / exceed585$total_area
exceed585$pct_exceed3_ph_area <- exceed585$ph_exceed3_area / exceed585$total_area
exceed585$pct_exceed_dual_area <- exceed585$dual_exceed_area / exceed585$total_area




### A

# percentage area of overall coral exceeded envelope vs percent area of species exceeded envelopes
temp_envelope <-  ggplot() +
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
  scale_y_continuous(limits = c(0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) +
  labs(x = "Year",
       y = "Area exceeding \ntemperature envelopes (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )





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
  scale_y_continuous(limits = c(0, 5.5),breaks=(seq(0, 5, 1))) +
  labs(x = "Year",
       y = "Summer temperature \ndifference (°C)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )




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
    underestimated_pct_dual = pct_exceed_dual_area - allcoral_pct_dual
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
    dual_avg_magnitude = mean(underestimated_pct_dual[underestimated_pct_dual > 0], na.rm = TRUE) * 100
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
                y = mean_underestimate*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area underestimated by \nhabitat-level temperature analysis (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))+
  scale_y_continuous(limits = c(-0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) 




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
       y = "Area falling below \npH envelopes (%)") +
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
                y = mean_summer_delta,
                group = scenario),
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
       y = "Summer pH difference") +
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
                y = mean_underestimate*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # Zero reference line (where species = allcoral)
  geom_hline(yintercept = 0, color = "red", linewidth = 1) +
  labs(x = "Year",
       y = "Area underestimated by \nhabitat-level pH analysis (%)") +
  theme_classic() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA))+
  scale_y_continuous(limits = c(-0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) 









## G

#both temperature and ph envelopes exceeded 
dual_envelope <-  ggplot() +
  #species lines
  geom_line(data = exceed585 %>% filter(species_id != "allcoral"),
            aes(x = year, 
                y = pct_exceed_dual_area*100,
                group = interaction(species_id, scenario)), 
            color = "#4292c6", alpha = 0.03, linewidth = 0.15) +
  # average of all species 
  geom_line(data = exceed585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_envelope = mean(pct_exceed_dual_area, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_envelope*100,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  # all coral habitat line
  geom_line(data = exceed585 %>% filter(species_id == "allcoral"),
            aes(x = year, 
                y = pct_exceed_dual_area*100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100), 
                     #labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 100, 20)) +  scale_color_viridis_d() +
  labs(x = "Year",
       y = "Area exceeding both \ntemperature and pH envelopes (%)") +
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





fig1 <- plot_grid(temp_envelope, temp_delta, temp_envelope_underestimate,
                  ph_envelope, ph_delta, ph_envelope_underestimate,
                  dual_envelope, NULL,dual_envelope_underestimate,
          labels = c("A", "B","C","D","E","F","G","","H"), ncol = 3,nrow=3)


ggsave(paste0(figures_folder, "figure_1.pdf"), fig1, 
       width = 10, height = 10, dpi = 300)  

ggsave(paste0(figures_folder, "figure_1.png"), fig1, 
       width = 10, height = 10, dpi = 350)

















### FIG 2

#see range_sim.R














### FIG 3

#simulation plot with real data

#temperature and pH profiles






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
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Mean summer \ntemperature (°C)",
       color = "Scenario") +
  scale_x_continuous(limits=c(-40,45))
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

ph_lat <- ggplot(combined_lat585, aes(x = latitude, y = ph, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Mean summer \npH",
       color = "Scenario") +
  scale_x_continuous(limits=c(-40,45))
  #theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#temp increase v range size
tempdelta_range <- ggplot(exceedsp585[exceedsp585$year==2100,], aes(x = total_area/exceedall585$total_area[1], y = mean_temp_delta)) +
  #geom_point(alpha = 0.6, color = "darkred") +
  geom_pointdensity()+
  scale_color_gradient(low='#eff3ff',high="#2171b5") + 
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Mean temperature \nincrease (°C)")+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2))+
  theme(legend.position="none")


# % loss vs range size
# Sort the data by range size
sorted_data <- exceedsp585[exceedsp585$year==2100,] %>%
  arrange(total_area/exceedall585$total_area[1])
# Calculate cumulative percentage of species
sorted_data$cum_species_pct <- seq_len(nrow(sorted_data)) / nrow(sorted_data)

tempexceed_range <- ggplot(exceed585[exceed585$year==2100,], aes(x = total_area/exceedall585$total_area[1], y = pct_exceed_temp_area)) +
  geom_pointdensity()+
  scale_color_gradient(low='#eff3ff',high="#2171b5") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  # Add the cumulative species line
  geom_line(data = sorted_data, aes(x = total_area/exceedall585$total_area[1], y = cum_species_pct), color = "blue", linewidth = 1) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Habitat Loss \n(% Area Above Today's Max Temp)")+
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(scale = 100))+
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100))+
  theme(legend.position="none")


#distribution plot - sizes of ranges
range_size <- ggplot(exceedsp585[exceedsp585$year==2100,],aes(x=total_area/exceedall585$total_area[1]))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Range as proportion \nof domain size", 
       y = "Count")+
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                   labels = scales::percent_format(scale = 100))+
  scale_y_continuous(limits=c(0,100))

#ph increase v range size
phdelta_range <- ggplot(exceedsp585[exceedsp585$year==2100,], aes(x = total_area/exceedall585$total_area[1], y = mean_ph_delta)) +
  #geom_point(alpha = 0.6, color = "darkred") +
  geom_pointdensity()+
  scale_color_gradient(low='#eff3ff',high="#2171b5") + 
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Mean ph decrease")+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2))+
  theme(legend.position="none")


# % loss vs range size
phexceed_range <- ggplot(exceed585[exceed585$year==2100,], aes(x = total_area/exceedall585$total_area[1], y = pct_exceed_ph_area)) +
  #geom_point(alpha = 0.6, color = "darkred") +
  geom_pointdensity()+
  scale_color_gradient(low='#eff3ff',high="#2171b5") + #scale_color for point, scale_fill for hex
  #scale_fill_viridis_c()+
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Habitat loss \n(% area below today's min pH)")+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2))+
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100))+
  theme(legend.position="none")

fig3_coral <- plot_grid(temp_lat, tempdelta_range, tempexceed_range,range_size, ph_lat,phdelta_range, phexceed_range,
          labels = c( "A","B","C","D","E","F","G"), ncol = 4,nrow=2)


ggsave(paste0(figures_folder, "figure_3.pdf"), fig3_coral, 
       width = 10, height = 10, dpi = 300)  

ggsave(paste0(figures_folder, "figure_3.png"), fig3_coral, 
       width = 16, height = 6, dpi = 350)



















### FIG 4
# range size v distance to analog

species_coord <- fread(paste0(output_directory, 'species_pixels_1982-1992.csv'))

diffanalogN_585 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeansummerN', "_", "ssp585", "_rangeTemp", 
                                                        1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))
diffanalogS_585 <- fread(paste0(output_directory_analog, "split/analogdiffdf_", 'modelmeansummerS', "_", "ssp585", "_rangeTemp", 
                                1, "_rangepH", 0.1, "_", 'hist', "-", 2100, ".csv"))


species_coord$px_temp <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
species_coord[, py_temp := y] 

north_result <- merge(species_coord[y >= 0], diffanalogN_585, by = c("px_temp", "py_temp"), all.x = TRUE)
south_result <- merge(species_coord[y < 0], diffanalogS_585, by = c("px_temp", "py_temp"), all.x = TRUE)
species_analog585 <- rbindlist(list(north_result, south_result), fill = TRUE)

sp_rangesize <- unique(exceedsp585[, .(species_id, total_area)])
sp_rangesize$species_id <- as.character(sp_rangesize$species_id)
species_analog585$species_id <- as.character(species_analog585$species_id)
species_analog585 <- merge(species_analog585, sp_rangesize, by = "species_id", all.x = TRUE)


temp_plot <- ggplot(species_analog585, aes(x = total_area, y = distance_temp)) +
  geom_point(size=0.1,alpha = 0.1) +  # Use alpha for better visibility if points overlap
  #scale_x_log10() +  # Log scale for total_area as range sizes often span orders of magnitude
  labs(x = "Range Size (km2)",
    y = "Distance to Temperature Analog") +
  theme_classic() 

ph_plot <- ggplot(species_analog585, aes(x = total_area, y = distance_ph)) +
  geom_point(size=0.1,alpha = 0.1) +
  #scale_x_log10() +  # Log scale for total_area
  labs(x = "Range Size (km2)",
    y = "Distance to pH Analog") +
  theme_classic()  

species_avg <- species_analog585[, .(
  avg_distance_temp = mean(distance_temp, na.rm = TRUE),
  se_distance_temp = sd(distance_temp, na.rm = TRUE) / sqrt(.N),
  avg_distance_ph = mean(distance_ph, na.rm = TRUE),
  se_distance_ph = sd(distance_ph, na.rm = TRUE) / sqrt(.N),
  total_area = first(total_area)  # Take the first total_area value for each species
), by = species_id]

temp_avg_plot <- ggplot(species_avg, aes(x = total_area, y = avg_distance_temp)) +
  geom_point(size = 0.1, alpha = 0.1) +  # Small points with high transparency
  #scale_x_log10() +  # Log scale for total_area
  labs(x = "Range Size (km2)",
    y = "Average Distance to Temperature Analog") +
  theme_classic() 

ph_avg_plot <- ggplot(species_avg, aes(x = total_area, y = avg_distance_ph)) +
  geom_point(size = 0.1, alpha = 0.1) +  # Small points with high transparency
  #scale_x_log10() +  # Log scale for total_area
  labs(x = "Range Size (km2)",
    y = "Average Distance to pH Analog") +
  theme_classic()

temp_error_plot <- ggplot(species_avg, aes(x = total_area, y = avg_distance_temp)) +
  geom_point(size = 0.1, alpha = 0.2) +  # Small points with high transparency
  geom_errorbar(aes(ymin = avg_distance_temp - se_distance_temp, 
                    ymax = avg_distance_temp + se_distance_temp),
                width = 0.1, alpha = 0.2) +  # Error bars with transparency
  #scale_x_log10() +  # Log scale for total_area
  labs(x = "Range Size (total area)",
    y = "Average Distance to Temperature Analog") +
  theme_classic() 

ph_error_plot <- ggplot(species_avg, aes(x = total_area, y = avg_distance_ph)) +
  geom_point(size = 0.1, alpha = 0.2) +  # Small points with high transparency
  geom_errorbar(aes(ymin = avg_distance_ph - se_distance_ph, 
                    ymax = avg_distance_ph + se_distance_ph),
                width = 0.1, alpha = 0.2) +  # Error bars with transparency
  #scale_x_log10() +  # Log scale for total_area
  labs(x = "Range Size (total area)",
    y = "Average Distance to pH Analog") +
  theme_classic() 


plot_grid(temp_plot, ph_plot,temp_avg_plot, ph_avg_plot,
          labels = c( "A","B","C","D"), ncol = 2,nrow=2)

plot_grid(temp_error_plot, ph_error_plot,
          labels = c( "A","B"), ncol = 2,nrow=1)












### FIG 5

#range size v habitat loss by 2100, colored by family

species_family <- species_envelope[, .(species_id = id_no, family)]
exceed585 <- merge(exceed585, species_family, by = "species_id")

ggplot(exceed585[exceed585$year==2100,], aes(x = total_area/exceedall585$total_area[1], y = pct_exceed_temp_area, color = family)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(x = "Range as proportion \nof domain size", 
       y = "Habitat Loss \n(% Area Above Today's Max Temp)")+
  scale_x_continuous(limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2),
                     labels = scales::percent_format(scale = 100))+
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100))+
  theme(legend.position="none")





#histogram by family
ggplot(exceedsp585[exceedsp585$year==2100,],aes(x=family))+
  geom_bar(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Family", 
       y = "Count")+
  scale_y_continuous(limits=c(0,60))+
  coord_flip()



# coral family and exceeded envelopes

species_envelope <- fread(paste0(output_directory, 'climate_envelopes_1982-1992.csv'))
species_family <- species_envelope[, .(species_id = id_no, family)]

all_species_data <- merge(all_species_data, species_family, by = "species_id")
exceed585 <- merge(exceed585,species_family,by = "species_id")

species_medians <- all_species_data[, .(
  temp_median = median(temp_exceed),
  ph_median = median(abs(ph_exceed)) 
), by = family]

temp_order <- species_medians[order(-temp_median)]$family  
ph_order <- species_medians[order(-ph_median)]$family

all_species_data[, temp_family := factor(family, levels = temp_order)]
all_species_data[, ph_family := factor(family, levels = ph_order)]

#as boxplot
fam_temp <- ggplot(all_species_data, aes(x = factor(family, levels = temp_order), 
                                         y = temp_exceed)) +  
  geom_boxplot(fill = "#FFA07A", alpha = 0.7) +
  geom_boxplot(data = all_coral_data,
               aes(x = "All Coral", y = temp_exceed),
               fill = "red", alpha = 0.7) +
  scale_x_discrete(limits = c("All Coral", temp_order)) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Family",
       y = "Temperature Exceedance (°C)") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

fam_ph <- ggplot(all_species_data, aes(x = factor(family, levels = ph_order), 
                                       y = ph_exceed)) +  
  geom_boxplot(fill = "#87CEEB", alpha = 0.7) +
  geom_boxplot(data = all_coral_data,
               aes(x = "All Coral", y = ph_exceed),
               fill = "blue", alpha = 0.7) +
  scale_x_discrete(limits = c("All Coral", ph_order)) +
  coord_flip() +
  theme_minimal() +
  labs(x = "Family",
       y = "pH Exceedance") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

plot_grid(fam_temp, fam_ph, labels = c("A", "B"), ncol = 2)















#as scatterplot
ggplot(all_species_data, aes(x = temp_exceed, y = ph_exceed, color = family)) +
  geom_point(alpha = 0.2, size = 2) +
  theme_minimal() +
  labs(x = "Temperature exceedance (°C)",
       y = "pH exceedance") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.title = element_blank())


#as scatterplot means only per species
species_means <- all_species_data[, .(
  mean_temp_exceed = mean(temp_exceed, na.rm = TRUE),
  mean_ph_exceed = mean(ph_exceed, na.rm = TRUE),
  total_sp_area = sum(sp_area,na.rm=T)
), by = .(species_id, family)]

ggplot(species_means, aes(x = mean_temp_exceed, y = mean_ph_exceed, color = total_sp_area)) +
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_viridis_c(option = "plasma", direction = -1) + 
  theme_minimal() +
  labs(x = "Mean temperature exceedance (°C)",
       y = "Mean pH exceedance",
       color = "Species range \nsize (km²)") +
  theme(panel.grid = element_blank(),
        legend.position = "right")

dhw585 <- dhw585 %>% mutate(species_id = as.character(species_id))
exceed585_2100 <- exceed585 %>% filter(year == 2100)
dhw_2100 <- dhw585 %>% filter(year == 2100)
combined_data <- exceed585_2100 %>%
  select(species_id, total_area, temp_exceed_area, ph_exceed1_area, dual_exceed_area,
         pct_exceed_temp_area, pct_exceed_ph_area, pct_exceed_dual_area) %>%
  left_join(dhw_2100 %>% select(species_id, dhw_exceed_area, dhw_exceed_pct),
            by = "species_id") %>%
  left_join(all_species_data %>% 
              select(species_id, family) %>% 
              distinct(),
            by = "species_id")



#as scatterplot means only per family
family_means <- exceed585_2100 %>%
  group_by(family) %>%
  summarize(
    mean_pct_temp_exceed_area = mean(pct_exceed_temp_area, na.rm = TRUE),
    mean_pct_ph_exceed_area = mean(pct_exceed_ph_area, na.rm = TRUE),
    mean_temp_exceed = mean(mean_temp_exceed, na.rm = TRUE),
    mean_ph_exceed = mean(mean_ph_exceed, na.rm = TRUE),
    mean_temp_diff = mean(mean_temp_diff, na.rm = TRUE),
    mean_ph_diff = mean(mean_ph_diff, na.rm = TRUE)
  )
ggplot() +
  # individual species points 
  geom_point(data = exceed585_2100, 
             aes(x = pct_exceed_temp_area, y = pct_exceed_ph_area, color = family),
             alpha = 0.2) +
  # family means 
  geom_point(data = family_means, 
             aes(x = mean_pct_temp_exceed_area, y = mean_pct_ph_exceed_area, color = family),
             size = 4) +
  #coral habitat
  geom_point(data=exceedall585[exceedall585$year==2100,],
             aes(x = temp_exceed_area/total_area, y = ph_exceed_area/total_area),
                 color="black", size=4)+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2)) +
  
  scale_y_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    x = "% area temperature exceeded",
    y = "% area pH exceeded",
    color = "Family"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#scatterplot exceedance temp & ph
ggplot() +
  # individual species points 
  geom_point(data = exceed585_2100, 
             aes(x = mean_temp_exceed, y = mean_ph_exceed, color = family),
             alpha = 0.2) +
  # family means 
  geom_point(data = family_means, 
             aes(x = mean_temp_exceed, y = mean_ph_exceed, color = family),
             size = 4) +
  #coral habitat
  geom_point(data=exceedall585[exceedall585$year==2100,],
             aes(x = mean_temp_exceed, y = mean_ph_exceed),
             color="black", size=4)+
  labs(
    x = "Mean temperature exceeded",
    y = "Mean pH exceeded",
    color = "Family"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#scatterplot difference temp & ph
ggplot() +
  # individual species points 
  geom_point(data = exceed585_2100, 
             aes(x = mean_temp_diff, y = mean_ph_diff, color = family),
             alpha = 0.2) +
  # family means 
  geom_point(data = family_means, 
             aes(x = mean_temp_diff, y = mean_ph_diff, color = family),
             size = 4) +
  #coral habitat
  geom_point(data=exceedall585[exceedall585$year==2100,],
             aes(x = mean_temp_diff, y = mean_ph_diff),
             color="black", size=4)+
  labs(
    x = "Mean temperature \n envelope vs value",
    y = "Mean pH \n envelope vs value",
    color = "Family"
  ) +
  theme_minimal() +
  theme(legend.position = "right")


#line dual stressor per year by family
family_dual <- exceed585 %>%
  group_by(family, year) %>%
  summarize(avg_dual_exceed_area = mean(dual_exceed_area, na.rm = TRUE),
            avg_pct_dual_exceed_area = mean(pct_exceed_dual_area,na.rm=T),
            .groups = "drop")

ggplot() +
  geom_line(data=family_dual,aes(x = year, y = avg_pct_dual_exceed_area, color = family, group = family),linewidth = 1) +
  geom_line(data=exceedall585,aes(x=year,y=dual_exceed_area/total_area),color="red",linewidth=2)+
  theme_minimal() +
  labs(x = "Year",
       y = "Mean % area \ndual exceedance",
       color = "Family") +
  theme(legend.position = "right",
        panel.grid.minor = element_blank(),
        text = element_text(size = 12))








# dhw asb year vs year
dhw585 <- rbind(dhwsp585,dhwall585)
dhw585$dhw_exceed_pct <- dhw585$dhw_exceed_area / dhw585$total_area
dhwsp585$dhw_exceed_pct <- dhwsp585$dhw_exceed_area / dhwsp585$total_area
dhwall585$dhw_exceed_pct <- dhwall585$dhw_exceed_area / dhwall585$total_area

dhw_asbarea <- ggplot() +
  geom_line(data=dhw585 %>% filter(species_id != "allcoral"),
            aes(x = year, y = dhw_exceed_pct, group = species_id),
            alpha = 0.1, linewidth = 0.2) + 
  # average of all species 
  geom_line(data = dhw585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_dhw_exceed = mean(dhw_exceed_pct, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_dhw_exceed,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  geom_line(data=dhw585%>% filter(species_id == "allcoral"),
            aes(x = year, y = dhw_exceed_pct),
            color="red",linewidth=1.5)+
  theme_minimal() +
  labs(x = "Year",
       y = "% area experiencing bleaching") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

dhw_asbdur <- ggplot() +
  geom_line(data=dhw585 %>% filter(species_id != "allcoral"),
            aes(x = year, y = mean_dhw_duration, group = species_id),
            alpha = 0.1, linewidth = 0.2) + 
  # average of all species 
  geom_line(data = dhw585 %>% 
              filter(species_id != "allcoral") %>%
              group_by(year, scenario) %>%
              summarize(mean_dhw_duration = mean(mean_dhw_duration, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_dhw_duration,
                group = scenario),
            color = "blue",
            linewidth = 1.5) +
  geom_line(data=dhw585%>% filter(species_id == "allcoral"),
            aes(x = year, y = mean_dhw_duration),
            color="red",linewidth=1.5)+
  theme_minimal() +
  labs(x = "Year",
       y = "Mean Duration (months)") +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

plot_grid(dhw_asbarea, dhw_asbdur, labels = c("A", "B"), ncol = 2,nrow=1)


bleachingsp50 <- dhwsp585 %>%
  group_by(species_id, family) %>%
  filter(dhw_exceed_pct >= 0.50) %>%
  arrange(year) %>%
  slice(1) %>%
  select(species_id, family, year)
all_bleaching50 <- dhwall585 %>%
  filter(dhw_exceed_pct >= 0.50) %>%
  arrange(year) %>%
  slice(1) %>%
  mutate(family = "HABITAT") %>%
  select(species_id, family, year)
bleaching50 <- bind_rows(bleachingsp50, all_bleaching50)
family_counts <- bleaching50 %>%
  group_by(family) %>%
  summarize(count = n())

ggplot(bleaching50, aes(x = reorder(family, year, FUN = median), y = year)) +
  geom_boxplot(fill = "coral3", alpha = 0.8) +
  geom_boxplot(data = bleaching50 %>% filter(family == "ALL CORAL"), 
               fill = "red", alpha = 0.8) +  # Highlight ALL CORAL
  geom_text(data = family_counts, aes(x = family, y = 2020, label = count), 
            hjust = -0.2, size = 3) +
  scale_y_continuous(limits = c(2020, 2065)) +
  coord_flip() +
  labs(x = "Family",
       y = "Year") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

### supplements











#quartiles 
#
species_sizes <- exceed585 %>%
  filter(species_id != "allcoral" & year == 2100) %>%
  select(species_id, total_area) %>%
  distinct()

# Calculate quartiles and assign each species to a quartile
quartile_breaks <- quantile(species_sizes$total_area, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
species_sizes <- species_sizes %>%
  mutate(size_quartile = cut(total_area, 
                             breaks = quartile_breaks,
                             labels = c("Q1 (Small)", "Q2 (Medium-Small)", "Q3 (Medium-Large)", "Q4 (Large)"),
                             include.lowest = TRUE))






# Join this information back to the exceed585_long dataset
exceed585_long_with_quartiles <- exceed585_long %>%
  left_join(species_sizes %>% select(species_id, size_quartile), by = "species_id")

temp_envelope_quart <-  ggplot() +
  #species lines
  geom_line(data = exceed585_long_with_quartiles %>% filter(species_id != "allcoral" & variable == "Temperature"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = interaction(species_id, scenario),
            ), #color = scenario
            alpha = 0.1, linewidth = 0.2) +
  # quartile mean lines
  geom_line(data = exceed585_long_with_quartiles %>% 
              filter(species_id != "allcoral" & variable == "Temperature" & !is.na(size_quartile)) %>%
              group_by(year, scenario, size_quartile) %>%
              summarize(mean_unsuitable = mean(unsuitable_area, na.rm = TRUE), .groups = "drop"),
            aes(x = year, 
                y = mean_unsuitable * 100,
                group = interaction(size_quartile, scenario),
                color = size_quartile),
            linewidth = 1.2) +
  # all coral habitat line
  geom_line(data = exceed585_long %>% filter(species_id == "allcoral" & variable == "Temperature"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = scenario),
            color = "red",
            linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("#ADD8E6", "#6495ED", "#4169E1", "#00008B"),  
                     name = "Range Size Quartile") +
  labs(x = "Year",
       y = "% area exceeding temperature envelopes") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )





legend_data <- data.frame(
  x = rep(1, 8),
  y = 1:8,
  variable = rep(c("pH", "pH1", "pH2", "pH3"), each = 2),
  coral_type = rep(c("All coral", "Individual species"), 4)
)

# Create color and linetype scales
ph_label_map <- c(
  "pH" = "envelope exceeded",
  "pH1" = "0.1 pH below",
  "pH2" = "0.2 pH below", 
  "pH3" = "0.3 pH below"
)

color_map <- c(
  "pH" = "black",
  "pH1" = "#FFEA00",
  "pH2" = "orange",
  "pH3" = "red"
)

ph_envelope_legend <-  ggplot() +
  #species lines
  geom_line(data = exceed585_long %>% filter(species_id != "allcoral" & variable == "pH"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = interaction(species_id, scenario)),
            color = "black",
            alpha = 0.1, linewidth = 0.2) +
  # all coral habitat line
  geom_line(data = exceed585_long %>% filter(species_id == "allcoral" & variable == "pH"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = scenario),
            color = "black",
            linewidth = 1.5,linetype="dashed") +
  #species lines
  geom_line(data = exceed585_long %>% filter(species_id != "allcoral" & variable == "pH1"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = interaction(species_id, scenario)),
            color = "#FFEA00", 
            alpha = 0.1, linewidth = 0.2) +
  # all coral habitat line
  geom_line(data = exceed585_long %>% filter(species_id == "allcoral" & variable == "pH1"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = scenario),
            color = "#FFEA00",
            linewidth = 1.5,linetype="dashed") +
  #species lines
  geom_line(data = exceed585_long %>% filter(species_id != "allcoral" & variable == "pH2"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = interaction(species_id, scenario)),
            color = "orange", 
            alpha = 0.1, linewidth = 0.2) +
  # all coral habitat line
  geom_line(data = exceed585_long %>% filter(species_id == "allcoral" & variable == "pH2"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = scenario),
            color = "orange",
            linewidth = 1.5,linetype="dashed") +
  #species lines
  geom_line(data = exceed585_long %>% filter(species_id != "allcoral" & variable == "pH3"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = interaction(species_id, scenario)),
            color = "red", 
            alpha = 0.1, linewidth = 0.2) +
  # all coral habitat line
  geom_line(data = exceed585_long %>% filter(species_id == "allcoral" & variable == "pH3"),
            aes(x = year, 
                y = unsuitable_area * 100,
                group = scenario),
            color = "red",
            linewidth = 1.5,linetype="dashed") +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_viridis_d() +
  labs(x = "Year",
       y = "Area exceeding pH envelopes (%)") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA)
  )+
  # Add a geom_line with dummy data solely for the legend
  geom_line(data = legend_data,
            aes(x = x, y = y, 
                color = variable,
                linetype = coral_type,
                linewidth = coral_type),
            alpha = 1,
            show.legend = TRUE) +
  scale_color_manual(values = color_map,
                     labels = ph_label_map,
                     name = "") +
  scale_linetype_manual(values = c("All coral" = "dashed", 
                                   "Individual species" = "solid"),
                        name = "") +
  scale_linewidth_manual(values = c("All coral" = 1.5, 
                                    "Individual species" = 0.8),
                         guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Year",
       y = "Area exceeding pH envelopes (%)") +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "right"
  ) +
  # Hide the dummy geom_line
  coord_cartesian(xlim = range(exceed585_long$year))

ph_legend <- get_legend(ph_envelope_legend)

ph_envelope <- ph_envelope_legend +
  theme(legend.position = "none")






