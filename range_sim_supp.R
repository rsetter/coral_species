


# variation in species frequency distributions

# Define all distribution types to compare
distribution_types <- c("uniform", "normal", "left_skewed", "right_skewed", "lognormal", "inverse_lognormal")

# Create a temperature gradient
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))
future_temp <- temperature_gradient(max_temp=33, min_temp=3, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()

# Run the simulation for each distribution type
for (dist_type in distribution_types) {
  segments <- generate_segments(1000, 1, 1000, latitude_range = c(0, 1000), distribution = dist_type)
  results <- analyze_temp_exceed(segments, today_temp, future_temp)
  
  # Add distribution type to results
  results$distribution <- dist_type
  
  results$pct_domain <- results$size / 1000
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
}

# Create a nicer label mapping for the plot legend
label_map <- c(
  "uniform" = "Uniform",
  "normal" = "Normal",
  "left_skewed" = "Left Skewed",
  "right_skewed" = "Right Skewed",
  "lognormal" = "Lognormal",
  "inverse_lognormal" = "Inverse Lognormal"
)
all_results$distribution_label <- label_map[all_results$distribution]

# Plot all trend lines using stat_summary_bin
comparison_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = distribution_label)) +
  #stat_summary_bin(fun = "mean", geom = "line",linewidth = 1.2,bins = 10) +
  geom_smooth(method = "loess", span = 0.3,se=F,linewidth=1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Habitat loss\n(% above today's max temp)",
       color = "Distribution Type") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 100))+
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = "none")

range_size_plot <- ggplot(all_results, aes(x = pct_domain*100, fill = distribution_label, color = distribution_label)) +
  geom_density(fill = NA, linewidth = 1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Density") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")  

legend_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = distribution_label)) +
  geom_line() +  # Use geom_line instead of geom_density for line symbols
  scale_color_brewer(palette = "Set1") +
  theme_classic() +
  labs(color = "Distribution type") +
  theme(legend.position = "right",  # Place legend on the right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        # Ensure legend uses lines not boxes
        legend.key = element_rect(fill = "white", color = NA))

legend <- cowplot::get_legend(legend_plot)

plot_grid(range_size_plot, comparison_plot, legend,labels = c("A", "B",""), ncol = 3,rel_widths = c(1, 1, 0.4))













## variation in temperature change

# Calculate temperature difference between Future and Today for each latitude
temp_changes585 <- combined_lat585 %>%
  select(latitude, temperature, scenario) %>%
  pivot_wider(names_from = scenario, values_from = temperature) %>%
  mutate(temp_change = Future - Today)
mean_tempdelta585 <- mean(temp_changes585$temp_change, na.rm = TRUE)
temp_changes370 <- combined_lat370 %>%
  select(latitude, temperature, scenario) %>%
  pivot_wider(names_from = scenario, values_from = temperature) %>%
  mutate(temp_change = Future - Today)
mean_tempdelta370 <- mean(temp_changes370$temp_change, na.rm = TRUE)
temp_changes245 <- combined_lat245 %>%
  select(latitude, temperature, scenario) %>%
  pivot_wider(names_from = scenario, values_from = temperature) %>%
  mutate(temp_change = Future - Today)
mean_tempdelta245 <- mean(temp_changes245$temp_change, na.rm = TRUE)
temp_changes126 <- combined_lat126 %>%
  select(latitude, temperature, scenario) %>%
  pivot_wider(names_from = scenario, values_from = temperature) %>%
  mutate(temp_change = Future - Today)
mean_tempdelta126 <- mean(temp_changes126$temp_change, na.rm = TRUE)

# Define the temperature increase scenarios to compare
temp_increases <- c(1, 3, 5, 10, 20,
                    mean_tempdelta126,mean_tempdelta245,mean_tempdelta370,mean_tempdelta585)

# Create the base temperature gradient for today
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()

# Run the simulation for each temperature increase scenario
for (temp_increase in temp_increases) {
  # Create future temperature with the specified increase
  future_temp <- temperature_gradient(
    max_temp = 30 + temp_increase, 
    min_temp = 0 + temp_increase,
    latitude_range = c(0, 1000)
  )
  
  # Generate segments using lognormal distribution
  segments <- generate_segments(
    1000, 1, 1000, 
    latitude_range = c(0, 1000), 
    distribution = "lognormal"
  )
  
  # Analyze the temperature exceedance
  results <- analyze_temp_exceed(segments, today_temp, future_temp)
  results$pct_domain <- results$size / 1000
  
  # Add temperature increase to results
  results$temp_increase <- temp_increase
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
}

all_results$scenario_label <- NA

# Set labels for the first 5 types (1-5)
all_results$scenario_label[all_results$temp_increase %in% temp_increases[1:5]] <- 
  paste0(temp_increases[1:5])[match(all_results$temp_increase[all_results$temp_increase %in% temp_increases[1:5]], 
                    temp_increases[1:5])]


# Set labels for the last 4 types (SSP scenarios)
all_results$scenario_label[all_results$temp_increase == temp_increases[6]] <- "ssp126"
all_results$scenario_label[all_results$temp_increase == temp_increases[7]] <- "ssp245"
all_results$scenario_label[all_results$temp_increase == temp_increases[8]] <- "ssp370"
all_results$scenario_label[all_results$temp_increase == temp_increases[9]] <- "ssp585"

# Convert to factor with specific order
all_results$scenario_label <- factor(
  all_results$scenario_label,
  levels = c("1", "3", "5", "10", "20", "ssp126", "ssp245", "ssp370", "ssp585")
) 

# Create a color palette that shows temperature intensity
temp_colors <- c("#99CCFF", "#3399FF", "#FF9933", "#FF6600", "#CC0000",
                 "lightgrey", "grey", "darkgrey", "black")



# Plot all trend lines using stat_summary_bin
comparison_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = scenario_label)) +
  #stat_summary_bin(fun = "mean", geom = "line",linewidth = 1.2,bins = 10) +
  geom_smooth(method = "loess", span = 0.3,se=F,linewidth=1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Habitat loss\n(% above today's max temp)",
       color = "Temperature\nIncrease") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_manual(values = temp_colors) +
  theme(
    legend.position = "right"
  )+
  theme(legend.position = "none")

temp_incr_plot <- ggplot(all_results, aes(x = pct_domain*100, y = temp_increase, color = scenario_label)) +
  geom_line(size = 1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Mean temperature \nincrease (°C)",
       color = "Temperature\nIncrease") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_y_continuous(limits = c(0, 20)) +  # Adjusted to show full range of temp increases (1-20)
  scale_color_manual(values = temp_colors) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")

legend_plot <- ggplot(all_results, aes(x = pct_domain*100, y = temp_increase, color = scenario_label)) +
  geom_line(size = 1) +
  scale_color_manual(values = temp_colors,
                     name = "Temperature\nIncrease") +
  theme_classic() +
  theme(legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(0.8, "cm"))

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

plot_grid(temp_incr_plot, comparison_plot, legend,labels = c("A", "B",""), ncol = 3,rel_widths = c(1, 1, 0.4))



















## variation in latitudinal temperature change



# Define the temperature change scenarios
scenarios <- list(
  "Uniform +3°C" = list(min_change = 3, max_change = 3),
  "Polar +5°C, Equator +0°C" = list(min_change = 5, max_change = 0),
  "Polar +5°C, Equator +3°C" = list(min_change = 5, max_change = 3),
  "Polar +0°C, Equator +5°C" = list(min_change = 0, max_change = 5),
  "Polar +3°C, Equator +5°C" = list(min_change = 3, max_change = 5)
)

# Create the base temperature gradient for today
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()

temp_increases <- data.frame()

# Run the simulation for each temperature change scenario
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  
  # Create future temperature with the specified changes
  future_temp <- temperature_gradient(
    max_temp = 30 + scenario$max_change, 
    min_temp = 0 + scenario$min_change,
    latitude_range = c(0, 1000)
  )
  
  increase_data <- data.frame(
    latitude = today_temp$latitude,
    temp_increase = future_temp$temperature - today_temp$temperature,
    scenario = scenario_name)
  temp_increases <- rbind(temp_increases, increase_data)
  
  
  # Generate segments using lognormal distribution
  segments <- generate_segments(
    1000, 1, 1000, 
    latitude_range = c(0, 1000), 
    distribution = "lognormal"
  )
  
  # Analyze the temperature exceedance
  results <- analyze_temp_exceed(segments, today_temp, future_temp)
  results$pct_domain <- results$size / 1000
  
  # Add scenario information to results
  results$scenario <- scenario_name
  results$min_change <- scenario$min_change
  results$max_change <- scenario$max_change
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
}


#  add temperature data scenario

# Create a smooth temperature function based on the actual data
create_temp_gradient_from_data <- function(lat_temp_data, latitude_range, smoothing = 0.7) {
  # Create a sequence of latitudes covering the full range
  lat_seq <- seq(latitude_range[1], latitude_range[2], length.out = 1000)
  
  # Calculate the midpoint in terms of the sequence index
  midpoint_index <- ceiling(length(lat_seq) / 2)
  midpoint_lat <- lat_seq[midpoint_index]
  
  # Normalize the actual latitude data to the range we're using
  min_actual_lat <- min(lat_temp_data$latitude)
  max_actual_lat <- max(lat_temp_data$latitude)
  actual_lat_range <- max_actual_lat - min_actual_lat
  
  # Use smooth.spline with adjustable smoothing parameter
  fit <- smooth.spline(lat_temp_data$latitude, lat_temp_data$temperature, spar = smoothing)
  
  # Create temperature array
  temps <- numeric(length(lat_seq))
  
  # Calculate temperatures for the first half (including midpoint)
  for (i in 1:midpoint_index) {
    # Convert simulation latitude to real-world latitude for lookup
    norm_pos <- (lat_seq[i] - latitude_range[1]) / (latitude_range[2] - latitude_range[1])
    real_lat <- min_actual_lat + norm_pos * actual_lat_range
    
    # Ensure we stay within the bounds of the real data
    real_lat <- max(min_actual_lat, min(max_actual_lat, real_lat))
    
    # Get temperature at this latitude using the smoothed fit
    temps[i] <- predict(fit, real_lat)$y
  }
  
  # Mirror the first half for the second half
  for (i in (midpoint_index+1):length(lat_seq)) {
    # Find the corresponding index in the first half
    mirror_index <- 2 * midpoint_index - i
    # Make sure we don't go below the first index
    mirror_index <- max(1, mirror_index)
    temps[i] <- temps[mirror_index]
  }
  
  # Create a temperature gradient object in the format expected by analyze_temp_exceed
  temp_gradient <- list(
    latitude = lat_seq,
    temperature = temps
  )
  
  return(temp_gradient)
}

today_data <- combined_lat585[combined_lat585$scenario == "Today",]
today_data <- today_data[order(today_data$latitude),]
today_temp <- create_temp_gradient_from_data(today_data, c(0, 1000))
today_df <- data.frame(latitude = today_temp$latitude,temperature = today_temp$temperature,scenario = "Today")

future_data <- combined_lat585[combined_lat585$scenario == "Future",]
future_data <- future_data[order(future_data$latitude),]
ssp585_future_temp <- create_temp_gradient_from_data(future_data, c(0, 1000))
ssp585_future_df <- data.frame(latitude = ssp585_future_temp$latitude, temperature = ssp585_future_temp$temperature,scenario = "Future")

future_data <- combined_lat370[combined_lat370$scenario == "Future",]
future_data <- future_data[order(future_data$latitude),]
ssp370_future_temp <- create_temp_gradient_from_data(future_data, c(0, 1000))
ssp370_future_df <- data.frame(latitude = ssp370_future_temp$latitude, temperature = ssp370_future_temp$temperature,scenario = "Future")

future_data <- combined_lat245[combined_lat245$scenario == "Future",]
future_data <- future_data[order(future_data$latitude),]
ssp245_future_temp <- create_temp_gradient_from_data(future_data, c(0, 1000))
ssp245_future_df <- data.frame(latitude = ssp245_future_temp$latitude, temperature = ssp245_future_temp$temperature,scenario = "Future")

future_data <- combined_lat126[combined_lat126$scenario == "Future",]
future_data <- future_data[order(future_data$latitude),]
ssp126_future_temp <- create_temp_gradient_from_data(future_data, c(0, 1000))
ssp126_future_df <- data.frame(latitude = ssp126_future_temp$latitude, temperature = ssp126_future_temp$temperature,scenario = "Future")


sim_clim585 <- rbind(today_df, ssp585_future_df)
sim_clim370 <- rbind(today_df, ssp370_future_df)
sim_clim245 <- rbind(today_df, ssp245_future_df)
sim_clim126 <- rbind(today_df, ssp126_future_df)


ggplot(sim_clim126, aes(x = latitude, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Mean summer \ntemperature (°C)",
       color = "Scenario") +
  scale_x_continuous(limits = c(0, 1000))

# Generate segments using lognormal distribution
segments <- generate_segments(1000, 1, 1000, latitude_range = c(0, 1000), distribution = "lognormal")

# Analyze temperature exceedance for actual data
ssp585_results <- analyze_temp_exceed(segments, today_temp, ssp585_future_temp)
ssp585_results$scenario <- "SSP585"
ssp585_results$min_change <- NA  
ssp585_results$max_change <- NA  
ssp585_results$pct_domain <- results$size / 1000

ssp370_results <- analyze_temp_exceed(segments, today_temp, ssp370_future_temp)
ssp370_results$scenario <- "SSP370"
ssp370_results$min_change <- NA  
ssp370_results$max_change <- NA  
ssp370_results$pct_domain <- results$size / 1000

ssp245_results <- analyze_temp_exceed(segments, today_temp, ssp245_future_temp)
ssp245_results$scenario <- "SSP245"
ssp245_results$min_change <- NA  
ssp245_results$max_change <- NA  
ssp245_results$pct_domain <- results$size / 1000

ssp126_results <- analyze_temp_exceed(segments, today_temp, ssp126_future_temp)
ssp126_results$scenario <- "SSP126"
ssp126_results$min_change <- NA  
ssp126_results$max_change <- NA  
ssp126_results$pct_domain <- results$size / 1000

# Combine with other results
all_results <- rbind(all_results, ssp585_results)
all_results <- rbind(all_results, ssp370_results)
all_results <- rbind(all_results, ssp245_results)
all_results <- rbind(all_results, ssp126_results)

#add ssp scenarios to tempincr
ssp585_increase <- data.frame(
  latitude = today_temp$latitude,
  temp_increase = ssp585_future_temp$temperature - today_temp$temperature,
  scenario = "SSP585"
)

ssp370_increase <- data.frame(
  latitude = today_temp$latitude,
  temp_increase = ssp370_future_temp$temperature - today_temp$temperature,
  scenario = "SSP370"
)

ssp245_increase <- data.frame(
  latitude = today_temp$latitude,
  temp_increase = ssp245_future_temp$temperature - today_temp$temperature,
  scenario = "SSP245"
)

ssp126_increase <- data.frame(
  latitude = today_temp$latitude,
  temp_increase = ssp126_future_temp$temperature - today_temp$temperature,
  scenario = "SSP126"
)

# Combine all temperature increase data
temp_increases <- rbind(
  temp_increases, 
  ssp585_increase, 
  ssp370_increase, 
  ssp245_increase, 
  ssp126_increase
)

temp_increases$latitude_degrees <- -50 + 100 * (temp_increases$latitude - min(temp_increases$latitude)) / 
  (max(temp_increases$latitude) - min(temp_increases$latitude))



# Set a custom color palette
scenario_colors <- c(
  "Uniform +3°C" = "#DC3912", "Polar +5°C, Equator +0°C" = "#FF9900", "Polar +5°C, Equator +3°C" = "#109618", 
  "Polar +0°C, Equator +5°C" = "#990099","Polar +3°C, Equator +5°C" = "#3366CC",
  "SSP126" = "lightgrey","SSP245" = "grey","SSP370" = "darkgrey","SSP585" = "black")

# Plot for habitat loss vs range size
comparison_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = scenario)) +
  #stat_summary_bin(fun = "mean", geom = "line",linewidth = 1.2,bins = 10) +
  geom_smooth(method = "loess", span = 0.3,se=F,linewidth=1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Habitat loss\n(% above today's max temp)",
       color = "Warming\nScenario") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) +
  scale_color_manual(values = scenario_colors) +
  theme(legend.position = "right")+
  theme(legend.position = "none")

# Plot for mean temperature increase vs range size
plot2 <-  ggplot(temp_increases, aes(x = latitude_degrees, y = temp_increase, color = scenario)) +
  geom_line(size = 1.2) +
  scale_color_manual(values = scenario_colors) +
  labs(
    x = "Latitude (degrees)",
    y = "Temperature increase (°C)",
    color = "Scenario"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_x_continuous(breaks = seq(-50, 50, by = 50))+
  theme(legend.position = "none")


legend_plot <- ggplot(temp_increases, aes(x = latitude_degrees, y = temp_increase, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors,
                     name = "Temperature\nIncrease") +
  theme_classic() +
  theme(legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(0.8, "cm"))

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

plot_grid(plot2, comparison_plot, legend,labels = c("A", "B",""), ncol = 3,rel_widths = c(1, 1, 0.4))

