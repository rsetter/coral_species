

#function for creating plots
create_plot_gridA <- function(today_temp, future_temp, results, analyzed_segments, displacement, 
                               labels = c("A", "B", "C", "D"), title = NULL, distribution_type = "uniform") {
  
  # Define color mapping based on distribution type
  color_map <- c(
    "uniform" = "#E41A1C",  # Red from Set1
    "normal" = "#377EB8",   # Blue from Set1
    "left_skewed" = "#4DAF4A",  # Green from Set1
    "right_skewed" = "#984EA3",  # Purple from Set1
    "lognormal" = "#FF7F00"  # Orange from Set1
  )
  
  # Combine temperature data
  combined_temp <- rbind(
    cbind(today_temp, scenario = "Today"),
    cbind(future_temp, scenario = "Future")
  )
  
  combined_temp <- combined_temp %>%
    mutate(display_latitude = ((latitude - min(latitude)) / (max(latitude) - min(latitude))) * 100 - 50)
  
  # Temperature by latitude plot
  plot1 <- ggplot(combined_temp, aes(x = display_latitude, y = temperature, color = scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
    theme_classic() +
    labs(x = "Latitude", 
         y = "Temperature (°C)",
         color = "Scenario") +
    scale_x_continuous(breaks = seq(-50, 50, by = 50))

  # Calculate percentage domain for the results
  results$pct_domain <- results$size / 1000
  
  # Plot range size vs temperature change
  plot2 <- ggplot(results, aes(x = pct_domain*100, y = temp_increase)) +
    geom_pointdensity(position = position_jitter(width = 0.01, height = 0.05), size = 0.01) +
    scale_color_gradient(low = '#cbecff', high = "#2171b5") + 
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Mean temperature \nincrease (°C)") +
    scale_x_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(0, 5)) +
    theme(legend.position = "none")
  
  # Calculate cumulative species percentage
  results <- results %>%
    arrange(pct_domain) %>%
    mutate(cum_species_pct = (row_number() / n()) * 100)
  
  # Plot range size vs habitat loss
  plot3 <- ggplot(results, aes(x = pct_domain*100, y = exceed_max_pct)) +
    geom_pointdensity(size=0.8) +
    scale_color_gradient(low = '#cbecff', high = "#2171b5", trans = "log10") +
    geom_smooth(method = "loess", span = 0.3, color = "black", se = F) +
    # Add the cumulative species line
    geom_line(aes(y = cum_species_pct * max(exceed_max_pct) / 100), 
              color = "blue", linewidth = 1) +
    scale_y_continuous(
      name = "Habitat loss",
      sec.axis = sec_axis(~. * 100 / max(results$exceed_max_pct), 
                          name = "Cumulative \nspecies (%)")
    ) +
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Habitat loss") +
    scale_x_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 20)) +
    theme(legend.position="none",
          axis.title.y.right = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "blue"))
  
  # Plot range size distribution with color from the distribution type
  plot4 <- ggplot(results, aes(x = pct_domain*100)) +
    geom_histogram(fill = color_map[distribution_type], color = "black") +
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Frequency \nof species") +
    scale_x_continuous(limits = c(-5, 100),
                       breaks = seq(0, 100, 20))
  
  # Combine plots with plot_grid
  combined_plot <- plot_grid(plot1, plot2, plot3, plot4, 
                             labels = labels, 
                             ncol = 4)
  
  # Add title if provided
  if (!is.null(title)) {
    title_grob <- ggdraw() + 
      draw_label(title, fontface = "bold", x = 0.1, y = 0.5)
    
    # Add the title above the plots
    combined_plot <- plot_grid(title_grob, combined_plot, 
                               ncol = 1, rel_heights = c(0.1, 1))
  }
  
  return(combined_plot)
}

create_plot_gridBC <- function(today_temp, future_temp, results, labels = c("A", "B", "C", "D"), title = NULL) {
  # Combine temperature data
  combined_temp <- rbind(
    cbind(today_temp, scenario = "Today"),
    cbind(future_temp, scenario = "Future")
  )
  
  combined_temp <- combined_temp %>%
    mutate(display_latitude = ((latitude - min(latitude)) / (max(latitude) - min(latitude))) * 100 - 50)
  
  # Temperature by latitude plot
  plot1 <- ggplot(combined_temp, aes(x = display_latitude, y = temperature, color = scenario)) +
    geom_line(size = 1) +
    scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
    theme_classic() +
    labs(x = "Latitude", 
         y = "Temperature (°C)",
         color = "Scenario") +
    scale_x_continuous(breaks = seq(-50, 50, by = 50))
  
  # Calculate percentage domain for the results
  results$pct_domain <- results$size / 1000
  
  # Plot range size vs temperature change
  plot2 <- ggplot(results, aes(x = pct_domain*100, y = temp_increase)) +
    geom_pointdensity(position = position_jitter(width = 0.01, height = 0.05), size = 0.01) +
    scale_color_gradient(low = '#cbecff', high = "#2171b5") + 
    geom_smooth(method = "loess", span = 0.3, color = "black", se = F,linewidth=0.5) +
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Mean temperature \nincrease (°C)") +
    scale_x_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 20)) +
    scale_y_continuous(limits = c(0, 5)) +
    theme(legend.position = "none")
  
  # Calculate cumulative species percentage
  results <- results %>%
    arrange(pct_domain) %>%
    mutate(cum_species_pct = (row_number() / n()) * 100)
  
  # Plot range size vs habitat loss
  plot3 <- ggplot(results, aes(x = pct_domain*100, y = exceed_max_pct)) +
    geom_pointdensity(size=0.8) +
    scale_color_gradient(low = '#cbecff', high = "#2171b5", trans = "log10") +
    geom_smooth(method = "loess", span = 0.3, color = "black", se = F) +
    # Add the cumulative species line
    geom_line(aes(y = cum_species_pct * max(exceed_max_pct) / 100), 
              color = "blue", linewidth = 1) +
    scale_y_continuous(
      name = "Habitat loss \n(% above today's max temp)",
      sec.axis = sec_axis(~. * 100 / max(results$exceed_max_pct), 
                          name = "Cumulative species (%)")
    ) +
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Habitat loss \n(% above today's max temp)") +
    scale_x_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, 20)) +
    theme(legend.position="none",
          axis.title.y.right = element_text(color = "blue"),
          axis.text.y.right = element_text(color = "blue"))
  
  # Plot range size distribution
  plot4 <- ggplot(results, aes(x = pct_domain*100)) +
    geom_histogram(fill = "lightgray", color = "black") +
    theme_classic() +
    labs(x = "Latitudinal range size (degrees)", 
         y = "Frequency \nof species") +
    scale_x_continuous(limits = c(-5, 100),
                       breaks = seq(0, 100, 20))
  
  # Combine plots with plot_grid
  combined_plot <- plot_grid(plot1, plot2, plot3, plot4, 
                             labels = labels, 
                             ncol = 4)
  
  # Add title if provided
  if (!is.null(title)) {
    title_grob <- ggdraw() + 
      draw_label(title, fontface = "bold", x = 0.1, y = 0.5)
    
    # Add the title above the plots
    combined_plot <- plot_grid(title_grob, combined_plot, 
                               ncol = 1, rel_heights = c(0.1, 1))
  }
  
  return(combined_plot)
}





# variation in species frequency distributions

# Define all distribution types to compare
distribution_types <- c("uniform", "normal", "left_skewed", "right_skewed", "lognormal")

# Create a temperature gradient
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))
future_temp <- temperature_gradient(max_temp=33, min_temp=3, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()
displacement_by_dist <- list()
analyzed_segments_by_dist <- list()

# Run the simulation for each distribution type
for (dist_type in distribution_types) {
  segments <- generate_segments(1000, 1, 1000, latitude_range = c(0, 1000), distribution = dist_type)
  results <- analyze_temp_exceed(segments, today_temp, future_temp)
  
  # Add distribution type to results
  results$distribution <- dist_type
  results$pct_domain <- results$size / 1000
  
  # Calculate latitudinal displacement
  displacement <- calculate_latitudinal_displacement(today_temp, future_temp)
  
  # Apply this to your segments
  analyzed_segments <- find_critical_range_size(displacement, segments)
  
  # Add distribution identifier to both datasets
  displacement$distribution <- dist_type
  analyzed_segments$distribution <- dist_type
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
  displacement_by_dist <- rbind(displacement_by_dist, displacement)
  analyzed_segments_by_dist <- rbind(analyzed_segments_by_dist, analyzed_segments)
}

# Create a nicer label mapping for the plot legend
label_map <- c(
  "uniform" = "Uniform",
  "normal" = "Normal",
  "left_skewed" = "Left Skewed",
  "right_skewed" = "Right Skewed",
  "lognormal" = "Lognormal"
)
library(RColorBrewer)
set1_colors <- brewer.pal(5, "Set1")
distribution_colors <- c(
  "uniform" = set1_colors[1],     # First color from Set1
  "normal" = set1_colors[2],      # Second color from Set1
  "left_skewed" = set1_colors[3], # Third color from Set1 
  "right_skewed" = set1_colors[4],# Fourth color from Set1
  "lognormal" = set1_colors[5]    # Fifth color from Set1
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
                     breaks = seq(0, 100, 20),
                     expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 100),
                     expand=c(0,0))+
  scale_color_manual(values = setNames(distribution_colors, label_map)) +
  theme(legend.position = "none")

range_size_plot <- ggplot(all_results, aes(x = pct_domain*100, fill = distribution_label, color = distribution_label)) +
  geom_density(fill = NA, linewidth = 1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Density") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.032),
                     expand=c(0,0.0005))+
  scale_fill_manual(values = setNames(distribution_colors, label_map)) +
  scale_color_manual(values = setNames(distribution_colors, label_map)) +
  theme(legend.position = "none")  


# Latitudinal displacement plot
lat_disp_distr <- ggplot(analyzed_segments, aes(x = display_center_lat, y = display_size)) +
  geom_line(data = displacement, aes(x = display_latitude, y = abs(display_displacement)), 
            color = "black", size = 1.5) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  scale_x_continuous(limits = c(-50, 50),
                     breaks = seq(-50, 50, 25),
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(0,100),
                     expand=c(0,0))

# percentage complete loss
complete_loss_dist <- analyzed_segments_by_dist %>%
  group_by(distribution) %>%
  summarize(
    total_count = n(),
    habitat_loss_count = sum(complete_habitat_loss),
    percentage = (habitat_loss_count / total_count) * 100)

complete_loss_dist_plot <-ggplot(complete_loss_dist, aes(x = distribution, y = percentage, fill = distribution)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = distribution_colors, labels = label_map) +  
  scale_x_discrete(labels = label_map) + 
  labs(x = "Distribution type",
    y = "Species with \ncomplete habitat loss (%)",
    fill = "Distribution Type") +
  scale_y_continuous(limits=c(0,50),expand=c(0,0))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none")

legend_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = distribution_label)) +
  geom_line(linewidth=1.2) +  # Use geom_line instead of geom_density for line symbols
  scale_color_manual(values = setNames(distribution_colors, label_map)) +
  theme_classic() +
  labs(color = "Distribution type") +
  theme(legend.position = "right",  # Place legend on the right
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        # Ensure legend uses lines not boxes
        legend.key = element_rect(fill = "white", color = NA))

legend <- cowplot::get_legend(legend_plot)

sim_distrib <- plot_grid(range_size_plot, comparison_plot, lat_disp_distr,complete_loss_dist_plot,legend,
                         labels = c("A", "B","C","D",""), ncol = 5,rel_widths = c(1, 1,1,1, 0.4),align="hv",axis = "tblr")

sim_uniform_dist <- create_plot_gridA(today_temp, future_temp,all_results[all_results$distribution =="uniform",],
                                      analyzed_segments_by_dist[analyzed_segments_by_dist$distribution =="uniform",],
                                      displacement_by_dist[displacement_by_dist$distribution == "uniform",],
                                      labels=c("E","F","G","H"), title="Uniform", distribution_type = "uniform")
sim_normal_dist <- create_plot_gridA(today_temp, future_temp,all_results[all_results$distribution =="normal",],
                                     analyzed_segments_by_dist[analyzed_segments_by_dist$distribution=="normal",],
                                     displacement_by_dist[displacement_by_dist$distribution=="normal",],
                                     labels=c("I","J","K","L"), title="Normal", distribution_type = "normal")
sim_left_dist <- create_plot_gridA(today_temp, future_temp,all_results[all_results$distribution =="left_skewed",],
                                   analyzed_segments_by_dist[analyzed_segments_by_dist$distribution=="left_skewed",],
                                   displacement_by_dist[displacement_by_dist$distribution =="left_skewed",],
                                   labels=c("M","N","O","P"), title="Left-skewed", distribution_type = "left_skewed")
sim_right_dist <- create_plot_gridA(today_temp, future_temp,all_results[all_results$distribution =="right_skewed",],
                                    analyzed_segments_by_dist[analyzed_segments_by_dist$distribution=="right_skewed",],
                                    displacement_by_dist[displacement_by_dist$distribution=="right_skewed",],
                                    labels=c("Q","R","S","T"), title="Right-skewed", distribution_type = "right_skewed")
sim_logn_dist <- create_plot_gridA(today_temp, future_temp,all_results[all_results$distribution =="lognormal",],
                                   analyzed_segments_by_dist[analyzed_segments_by_dist$distribution=="lognormal",],
                                   displacement_by_dist[displacement_by_dist$distribution=="lognormal",],
                                   labels=c("U","V","W","X"), title="Lognormal", distribution_type = "lognormal")

fig_s_distr <- plot_grid(
  sim_distrib, sim_uniform_dist, sim_normal_dist, sim_left_dist, sim_right_dist, sim_logn_dist,
  ncol = 1, rel_heights = c(1.3, 1, 1, 1, 1, 1))


ggsave(paste0(figures_folder, "fig_S_distribution.png"), fig_s_distr, 
       width = 12, height = 15, dpi = 350)





















## variation in temperature change

# Calculate temperature difference between Future and Today for each latitude
# Define the temperature increase scenarios to compare
temp_increases <- c(1, 2, 3, 4, 5)

# Create the base temperature gradient for today
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()
displacement_by_temp <- data.frame()
analyzed_segments_by_temp <- data.frame()

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
  
  # Add distribution type to results
  results$distribution <- dist_type
  results$pct_domain <- results$size / 1000
  
  # Calculate latitudinal displacement
  displacement <- calculate_latitudinal_displacement(today_temp, future_temp)
  
  # Apply this to your segments
  analyzed_segments <- find_critical_range_size(displacement, segments)
  
  # Add distribution identifier to both datasets
  displacement$temp_increase <- temp_increase
  analyzed_segments$temp_increase <- temp_increase
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
  displacement_by_temp <- rbind(displacement_by_temp, displacement)
  analyzed_segments_by_temp <- rbind(analyzed_segments_by_temp, analyzed_segments)
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
}

all_results$scenario_label <- NA

# Set labels for the first 5 types (1-5)
all_results$scenario_label[all_results$temp_increase %in% temp_increases[1:5]] <- 
  paste0(temp_increases[1:5])[match(all_results$temp_increase[all_results$temp_increase %in% temp_increases[1:5]], 
                    temp_increases[1:5])]



# Convert to factor with specific order
all_results$scenario_label <- factor(
  all_results$scenario_label,
  levels = c("1", "2", "3", "4", "5")
) 

# Create a color palette that shows temperature intensity
temp_colors <- c("#99CCFF", "#3399FF", "#FF9933", "#FF6600", "#CC0000")



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
  scale_y_continuous(limits = c(0, 5)) +  # Adjusted to show full range of temp increases (1-20)
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

# Plot Latitudinal displacement plot
lat_disp_temp <- ggplot(analyzed_segments_by_temp, aes(x = display_center_lat, y = display_size)) +
  # Regular points
  geom_point(color = "black", size = 0.5, alpha = 0.1) +
  # No analog points 
  geom_point(data = subset(analyzed_segments_by_temp, no_analog), 
             aes(x = display_center_lat, y = display_size),
             color = "black", size = 0.5, alpha = 0.1) +
  # required latitudinal shift
  geom_line(data = displacement_by_temp, 
            aes(x = display_latitude, y = abs(display_displacement), 
                color = factor(temp_increase), group = factor(temp_increase)), 
            size = 1.2) +
  scale_color_manual(values = c("#99CCFF", "#3399FF", "#FF9933", "#FF6600", "#CC0000"), 
                     labels = c("1°C", "3°C", "5°C", "10°C", "20°C"),
                     name = "Temperature \nIncrease") +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  theme(legend.position = "none")

# percentage complete loss
complete_loss_temp <- analyzed_segments_by_temp %>%
  group_by(temp_increase) %>%
  summarize(
    total_count = n(),
    habitat_loss_count = sum(complete_habitat_loss),
    percentage = (habitat_loss_count / total_count) * 100)
complete_loss_temp$temp_increase <- as.character(complete_loss_temp$temp_increase)

complete_loss_temp_plot <-ggplot(complete_loss_temp, aes(x = temp_increase, y = percentage, fill = temp_increase)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = temp_colors, labels = label_map) +  
  scale_x_discrete(labels = label_map) + 
  labs(x = "Temperature increase (°C)",
       y = "Species with \ncomplete habitat loss (%)",
       fill = "Temperature Increase") +
  scale_y_continuous(limits=c(0,60),expand=c(0,0))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

sim_temp <- plot_grid(temp_incr_plot, comparison_plot,lat_disp_temp, complete_loss_temp_plot,legend,
                      labels = c("A", "B","C","D",""), ncol = 5,rel_widths = c(1, 1,1, 1,0.4),align="hv",axis = "tblr")


sim_temp1 <- create_plot_gridBC(today_temp, future_temp,all_results[all_results$temp_increase ==1,],
                                      labels=c("E","F","G","H"), title="1°C")
sim_temp2 <- create_plot_gridBC(today_temp, future_temp,all_results[all_results$temp_increase ==2,],
                               labels=c("I","J","K","L"), title="2°C")
sim_temp3 <- create_plot_gridBC(today_temp, future_temp,all_results[all_results$temp_increase ==3,],
                                labels=c("M","N","O","P"), title="3°C")
sim_temp4 <- create_plot_gridBC(today_temp, future_temp,all_results[all_results$temp_increase ==4,],
                                    labels=c("Q","R","S","T"), title="4°C")
sim_temp5 <- create_plot_gridBC(today_temp, future_temp,all_results[all_results$temp_increase ==5,],
                                   labels=c("U","V","W","X"), title="5°C")

fig_s_temp <- plot_grid(
  sim_temp, sim_temp1, sim_temp2, sim_temp3, sim_temp4, sim_temp5,
  ncol = 1, rel_heights = c(1, 1, 1, 1, 1, 1))


ggsave(paste0(figures_folder, "fig_S_temp.png"), fig_s_temp, 
       width = 12, height = 15, dpi = 350)















## variation in latitudinal temperature change



# Define the temperature change scenarios
scenarios <- list(
  "Uniform +1°C" = list(min_change = 1, max_change = 1),
  "Polar +1°C, Equator +3°C" = list(min_change = 1, max_change = 3),
  "Polar +1°C, Equator +5°C" = list(min_change = 1, max_change = 5),
  "Polar +3°C, Equator +1°C" = list(min_change = 3, max_change = 1),
  "Polar +5°C, Equator +1°C" = list(min_change = 5, max_change = 1)
)



# Create the base temperature gradient for today
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()
displacement_by_shape <- data.frame()
analyzed_segments_by_shape <- data.frame()
temp_increases <- data.frame()

# Run the simulation for each temperature change scenario
for (scenario_name in names(scenarios)) {
  scenario <- scenarios[[scenario_name]]
  
  # Special case for the plateau scenario
  if (scenario_name == "Plateau +3°C") {
    # Create plateau temperature profiles
    today_temp_plateau <- temperature_gradient_plateau(
      max_temp = 30, 
      min_temp = 0,
      latitude_range = c(0, 1000), 
      plateau_range = c(-15, 15)
    )
    
    future_temp_plateau <- temperature_gradient_plateau(
      max_temp = 33, 
      min_temp = 3,
      latitude_range = c(0, 1000), 
      plateau_range = c(-15, 15)
    )
    
    # Use these special temperature profiles for this scenario only
    temp_today <- today_temp_plateau
    temp_future <- future_temp_plateau
  } else {
    # For all other scenarios, use the standard temperature gradient
    temp_today <- today_temp
    temp_future <- temperature_gradient(
      max_temp = 30 + scenario$max_change, 
      min_temp = 0 + scenario$min_change,
      latitude_range = c(0, 1000)
    )
  }
  
  # Calculate temperature increase for plotting
  increase_data <- data.frame(
    latitude = temp_today$latitude,
    temp_increase = temp_future$temperature - temp_today$temperature,
    scenario = scenario_name
  )
  temp_increases <- rbind(temp_increases, increase_data)
  
  # Rest of your code remains the same, using temp_today and temp_future
  # instead of directly using today_temp and future_temp
  
  # Generate segments using lognormal distribution
  segments <- generate_segments(
    1000, 1, 1000, 
    latitude_range = c(0, 1000), 
    distribution = "lognormal"
  )
  
  # Analyze the temperature exceedance
  results <- analyze_temp_exceed(segments, temp_today, temp_future)
  results$pct_domain <- results$size / 1000
  
  # Add scenario information to results
  results$scenario <- scenario_name
  results$min_change <- scenario$min_change
  results$max_change <- scenario$max_change
  
  # Calculate latitudinal displacement
  displacement <- calculate_latitudinal_displacement(temp_today, temp_future)
  
  # Apply this to your segments
  analyzed_segments <- find_critical_range_size(displacement, segments)
  
  # Add distribution identifier to both datasets
  displacement$scenario <- scenario_name
  analyzed_segments$scenario <- scenario_name
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
  displacement_by_shape <- rbind(displacement_by_shape, displacement)
  analyzed_segments_by_shape <- rbind(analyzed_segments_by_shape, analyzed_segments)
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
}


temp_increases$latitude_degrees <- -50 + 100 * (temp_increases$latitude - min(temp_increases$latitude)) / 
  (max(temp_increases$latitude) - min(temp_increases$latitude))

# Set a custom color palette
scenario_colors <- c(
  "Uniform +1°C" = "#109618", 
  "Polar +1°C, Equator +3°C" = "#FF9933","Polar +1°C, Equator +5°C" = "#CC0000",
  "Polar +3°C, Equator +1°C" = "#99CCFF", "Polar +5°C, Equator +1°C" = "#3399FF")

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
  geom_line(size = 1.2,alpha=0.8,position = position_nudge(y = ifelse(temp_increases$scenario == "Plateau +3°C", 0.02, 
                                                            ifelse(temp_increases$scenario == "Uniform +3°C", -0.02, 0)))) +
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
  scale_y_continuous(breaks = seq(0, 5, by = 1),
                     limits=c(0,5))+
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


legend_plot <- ggplot(temp_increases, aes(x = latitude_degrees, y = temp_increase, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors,
                     name = "Temperature\nIncrease",
                     labels = c("Polar +1°C, \nEquator +3°C", 
                                "Polar +1°C, \nEquator +5°C", 
                                "Polar +3°C, \nEquator +1°C", 
                                "Polar +5°C, \nEquator +1°C",
                                "Uniform +1°C")) +
  theme_classic() +
  theme(legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),  # Reduced text size
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(0.8, "cm"),  # Reduced key size
        legend.key.height = unit(0.8, "cm"),  # Reduced key height
        legend.spacing.y = unit(1, "cm"))   # Reduced spacing between legend items

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

# Plot Latitudinal displacement plot
lat_disp_shape <- ggplot(analyzed_segments_by_shape, aes(x = display_center_lat, y = display_size)) +
  # Regular points
  geom_point(color = "black", size = 0.5, alpha = 0.1) +
  # No analog points 
  geom_point(data = subset(analyzed_segments_by_shape, no_analog), 
             aes(x = display_center_lat, y = display_size),
             color = "black", size = 0.5, alpha = 0.1) +
  # required latitudinal shift
  geom_line(data = displacement_by_shape, 
            aes(x = display_latitude, y = abs(display_displacement), 
                color = factor(scenario), group = factor(scenario)), 
            size = 1.2) +
  scale_color_manual(values = c( "#FF9933","#CC0000", "#99CCFF","#3399FF",  "#109618")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  theme(legend.position = "none")

# percentage complete loss
complete_loss_shape <- analyzed_segments_by_shape %>%
  group_by(scenario) %>%
  summarize(
    total_count = n(),
    habitat_loss_count = sum(complete_habitat_loss),
    percentage = (habitat_loss_count / total_count) * 100)

complete_loss_shape_plot <-ggplot(complete_loss_shape, aes(x = scenario, y = percentage, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = scenario_colors, labels = label_map) +  
  scale_x_discrete(labels = label_map) + 
  labs(x = "Temperature \nincrease",
       y = "Species with \ncomplete habitat loss (%)",
       fill = "Spatial temperature Increase") +
  scale_y_continuous(limits=c(0,50),expand=c(0,0))+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")


sim_shape <- plot_grid(plot2, comparison_plot, lat_disp_shape,complete_loss_shape_plot,legend,
                       labels = c("A", "B","C","D",""), ncol = 5,rel_widths = c(1, 1, 1,1,0.5),align="hv",axis = "tblr")


future_temp1 <- temperature_gradient(max_temp = 31, min_temp = 1,latitude_range = c(0, 1000))
future_temp2 <- temperature_gradient(max_temp = 33, min_temp = 1,latitude_range = c(0, 1000))
future_temp3 <- temperature_gradient(max_temp = 35, min_temp = 1,latitude_range = c(0, 1000))
future_temp4 <- temperature_gradient(max_temp = 31, min_temp = 3,latitude_range = c(0, 1000))
future_temp5 <- temperature_gradient(max_temp = 31, min_temp = 5,latitude_range = c(0, 1000))

sim_shape1 <- create_plot_gridBC(today_temp, future_temp1,all_results[all_results$scenario == "Uniform +1°C",],
                                labels=c("E","F","G","H"), title="Uniform +1°C")
sim_shape2 <- create_plot_gridBC(today_temp, future_temp2,all_results[all_results$scenario =="Polar +1°C, Equator +3°C",],
                                labels=c("I","J","K","L"), title="Polar +1°C, Equator +3°C")
sim_shape3 <- create_plot_gridBC(today_temp, future_temp3,all_results[all_results$scenario =="Polar +1°C, Equator +5°C",],
                                labels=c("M","N","O","P"), title="Polar +1°C, Equator +5°C")
sim_shape4 <- create_plot_gridBC(today_temp, future_temp4,all_results[all_results$scenario =="Polar +3°C, Equator +1°C",],
                                 labels=c("Q","R","S","T"), title="Polar +3°C, Equator +1°C°C")
sim_shape5 <- create_plot_gridBC(today_temp, future_temp5,all_results[all_results$scenario =="Polar +5°C, Equator +1°C",],
                                 labels=c("U","V","W","X"), title="Polar +5°C, Equator +1°C")

fig_s_shape <- plot_grid(
  sim_shape, sim_shape1, sim_shape2, sim_shape3, sim_shape4, sim_shape5,
  ncol = 1, rel_heights = c(1, 1, 1, 1, 1, 1))


ggsave(paste0(figures_folder, "fig_S_shape.png"), fig_s_shape, 
       width = 14, height = 15, dpi = 350)




















# variation in species range positions

# Define all locations to compare
location_types <- c("antarctic", "capricorn", "equator", "cancer", "arctic")

location_bounds <- list(
  "antarctic" = c(0, 200),
  "capricorn" = c(200, 400),
  "equator" = c(400, 600),
  "cancer" = c(600, 800),
  "arctic" = c(800, 1000))

# Create a temperature gradient
today_temp <- temperature_gradient(max_temp=30, min_temp=0, latitude_range = c(0, 1000))
future_temp <- temperature_gradient(max_temp=33, min_temp=3, latitude_range = c(0, 1000))

# Initialize a dataframe to store all results
all_results <- data.frame()
displacement_by_loc <- list()
analyzed_segments_by_loc <- list()

# Run the simulation for each distribution type
for (loc_type in location_types) {
  # Get the bounds for this location
  bounds <- location_bounds[[loc_type]]
  
  # Generate segments within the specific latitude range for this location
  segments <- generate_segments(1000, 1, 200, latitude_range = bounds, distribution = "uniform")
  
  # Analyze temperature exceedance
  results <- analyze_temp_exceed(segments, today_temp, future_temp)
  
  # Add location type to results
  results$location <- loc_type
  results$pct_domain <- results$size / 1000  # Keep using full domain percentage for comparison
  
  # Calculate latitudinal displacement
  displacement <- calculate_latitudinal_displacement(today_temp, future_temp)
  
  # Apply this to your segments
  analyzed_segments <- find_critical_range_size(displacement, segments)
  
  # Add location identifier to both datasets
  displacement$location <- loc_type
  analyzed_segments$location <- loc_type
  
  # Combine with previous results
  all_results <- rbind(all_results, results)
  displacement_by_loc <- rbind(displacement_by_loc, displacement)
  analyzed_segments_by_loc <- rbind(analyzed_segments_by_loc, analyzed_segments)
}

# Create a nicer label mapping for the plot legend
label_map <- c(
  "antarctic" = "Antarctic",
  "capricorn" = "Capricorn",
  "equator" = "Equator",
  "cancer" = "Cancer",
  "arctic" = "Arctic")

location_colors <- c(
  "antarctic" = "#99CCFF",  
  "capricorn" = "#FF6600",  
  "equator" = "#CC0000",    
  "cancer" = "#FF9933",      
  "arctic" = "#3399FF"       
)
all_results$location_label <- label_map[all_results$location]

# Plot all trend lines 
comparison_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = location_label)) +
  geom_smooth(method = "loess", span = 0.3, se=F, linewidth=1.2,alpha=0.8,position = position_dodge(width = 3)) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Habitat loss\n(% above today's max temp)",
       color = "Location Type") +
  scale_x_continuous(limits = c(0, 101),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0)) +
  scale_y_continuous(limits = c(0, 101),
                     expand=c(0,0))+
  scale_color_manual(values = setNames(location_colors, label_map)) +
  theme(legend.position = "none")

lat_dens_plot <- ggplot(analyzed_segments_by_loc, aes(x = display_center_lat,color = location)) +
  geom_density(alpha = 0.3, linewidth = 1.2) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Density") +
  scale_x_continuous(limits = c(-50, 50),
                     breaks = seq(-50, 50, 25),
                     expand = c(0,0)) +
    scale_color_manual(values =c(
      "antarctic" = "#99CCFF",  
      "capricorn" = "#FF6600",  
      "equator" = "#CC0000",    
      "cancer" = "#FF9933",      
      "arctic" = "#3399FF" )) +
  theme(legend.position = "none")

latitude_breaks <- c(-50, -30, -10, 10, 30, 50)
lat_hist_plot <- ggplot(analyzed_segments_by_loc, aes(x = display_center_lat, fill = location)) +
  geom_histogram(breaks = latitude_breaks, alpha = 0.8) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Count") +
  scale_x_continuous(limits = c(-50, 50),
                     breaks = seq(-50, 50, 20),
                     expand = c(0,0)) +
  scale_fill_manual(values = c(
    "antarctic" = "#99CCFF",  
    "capricorn" = "#FF6600",  
    "equator" = "#CC0000",    
    "cancer" = "#FF9933",      
    "arctic" = "#3399FF" )) +
  theme(legend.position = "none")



# Latitudinal displacement plot
lat_disp_loc <- ggplot() +
  # Add points for segments colored by location
  geom_point(data = analyzed_segments_by_loc, 
             aes(x = display_center_lat, y = display_size, color = location),
             size = 0.5, alpha = 0.1) +
  # Add the displacement line on top
  geom_line(data = displacement_by_loc[displacement_by_loc$location == location_types[1],], 
            aes(x = display_latitude, y = abs(display_displacement)),
            color = "black", size = 1.5) +
  # Set colors for the locations
  scale_color_manual(values = c(
    "antarctic" = "#99CCFF",  
    "capricorn" = "#FF6600",  
    "equator" = "#CC0000",    
    "cancer" = "#FF9933",      
    "arctic" = "#3399FF" )) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  scale_x_continuous(limits = c(-50, 50),
                     breaks = seq(-50, 50, 25),
                     expand=c(0,0)) +
  scale_y_continuous(limits=c(0,100),
                     expand=c(0,0)) +
  theme(legend.position = "none")

# percentage complete loss
complete_loss_loc <- analyzed_segments_by_loc %>%
  group_by(location) %>%
  summarize(
    total_count = n(),
    habitat_loss_count = sum(complete_habitat_loss),
    percentage = (habitat_loss_count / total_count) * 100)

complete_loss_loc_plot <- ggplot(complete_loss_loc, aes(x = location, y = percentage, fill = location)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = location_colors, labels = label_map) +  
  scale_x_discrete(labels = label_map) + 
  labs(x = "Location type",
       y = "Species with \ncomplete habitat loss (%)",
       fill = "Location Type") +
  scale_y_continuous(limits=c(0,100), expand=c(0,0))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

legend_plot <- ggplot(all_results, aes(x = pct_domain*100, y = exceed_max_pct, color = location_label)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values = setNames(location_colors, label_map)) +
  theme_classic() +
  labs(color = "Location type") +
  theme(legend.position = "right",
        legend.box.margin = margin(0, 0, 0, 0),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.key = element_rect(fill = "white", color = NA))

legend <- cowplot::get_legend(legend_plot)

sim_location <- plot_grid(lat_hist_plot, comparison_plot, lat_disp_loc, complete_loss_loc_plot, legend,
                          labels = c("A", "B", "C", "D", ""), ncol = 5, rel_widths = c(1, 1, 1, 1, 0.4),
                          align="hv", axis = "tblr")

sim_antarctic <- create_plot_gridBC(today_temp, future_temp, 
                                   all_results[all_results$location == "antarctic",],
                                   labels=c("E", "F", "G", "H"), title="Antarctic")

sim_capricorn <- create_plot_gridBC(today_temp, future_temp,
                                   all_results[all_results$location == "capricorn",],
                                   labels=c("I", "J", "K", "L"), title="Capricorn")

sim_equator <- create_plot_gridBC(today_temp, future_temp,
                                 all_results[all_results$location == "equator",],
                                 labels=c("M", "N", "O", "P"), title="Equator")

sim_cancer <- create_plot_gridBC(today_temp, future_temp,
                                all_results[all_results$location == "cancer",],
                                labels=c("Q", "R", "S", "T"), title="Cancer")

sim_arctic <- create_plot_gridBC(today_temp, future_temp,
                                all_results[all_results$location == "arctic",],
                                labels=c("U", "V", "W", "X"), title="Arctic")

# Combine all into final figure
fig_s_location <- plot_grid(
  sim_location, sim_antarctic, sim_capricorn, sim_equator, sim_cancer, sim_arctic,
  ncol = 1, rel_heights = c(1, 1, 1, 1, 1, 1))

# Save the figure
ggsave(paste0(figures_folder, "fig_S_location.png"), fig_s_location, 
       width = 12, height = 15, dpi = 350)

























# summary plot showing displacement (distance to analog) and complete loss

loss_displace <- analyzed_segments_by_dist %>%
  # Filter out rows with NA displacement values
  filter(!is.na(display_displacement)) %>%
  # Group by distribution type
  group_by(distribution) %>%
  summarize(
    avg_displacement = mean(display_displacement, na.rm = TRUE),
    min_displacement = min(display_displacement, na.rm = TRUE),
    max_displacement = max(display_displacement, na.rm = TRUE),
    # Calculate percentage of species with complete habitat loss
    percent_habitat_loss = sum(complete_habitat_loss == TRUE, na.rm = TRUE) / n() * 100)

loss_temp <- analyzed_segments_by_temp %>%
  # Filter out rows with NA displacement values
  filter(!is.na(display_displacement)) %>%
  # Group by temperature increase
  group_by(temp_increase) %>%
  summarize(
    avg_displacement = mean(display_displacement, na.rm = TRUE),
    min_displacement = min(display_displacement, na.rm = TRUE),
    max_displacement = max(display_displacement, na.rm = TRUE),
    # Calculate percentage of species with complete habitat loss
    percent_habitat_loss = sum(complete_habitat_loss == TRUE, na.rm = TRUE) / n() * 100)
loss_temp$temp_increase <- as.character(loss_temp$temp_increase)

loss_shape <- analyzed_segments_by_shape %>%
  # Filter out rows with NA displacement values
  filter(!is.na(display_displacement)) %>%
  # Group by temperature increase
  group_by(scenario) %>%
  summarize(
    avg_displacement = mean(display_displacement, na.rm = TRUE),
    min_displacement = min(display_displacement, na.rm = TRUE),
    max_displacement = max(display_displacement, na.rm = TRUE),
    # Calculate percentage of species with complete habitat loss
    percent_habitat_loss = sum(complete_habitat_loss == TRUE, na.rm = TRUE) / n() * 100)

loss_loc <- analyzed_segments_by_loc %>%
  # Filter out rows with NA displacement values
  filter(!is.na(display_displacement)) %>%
  # Group by location factor
  group_by(location) %>%
  summarize(
    avg_displacement = mean(display_displacement, na.rm = TRUE),
    min_displacement = min(display_displacement, na.rm = TRUE),
    max_displacement = max(display_displacement, na.rm = TRUE),
    # Calculate percentage of species with complete habitat loss
    percent_habitat_loss = sum(complete_habitat_loss == TRUE, na.rm = TRUE) / n() * 100)


# Create the plot with single colors per group but different colors between groups
displace_loss <- ggplot() +
  # Distribution group
  geom_segment(data = loss_displace, 
               aes(x = min_displacement, xend = max_displacement, 
                   y = percent_habitat_loss, yend = percent_habitat_loss,
                   group = distribution),
               color = "#CC0000", 
               size = 0.5, alpha = 0.8) +
  geom_point(data = loss_displace,
             aes(x = avg_displacement, y = percent_habitat_loss,
                 group = distribution),
             color = "#CC0000",
             size = 2, alpha = 0.8) +
  
  # Temperature group 
  geom_segment(data = loss_temp, 
               aes(x = min_displacement, xend = max_displacement, 
                   y = percent_habitat_loss, yend = percent_habitat_loss,
                   group = temp_increase),
               color = "#FF9933", 
               size = 0.5, alpha = 0.8) +
  geom_point(data = loss_temp,
             aes(x = avg_displacement, y = percent_habitat_loss,
                 group = temp_increase),
             color = "#FF9933", 
             size = 2, alpha = 0.8) +
  
  # Scenario group 
  geom_segment(data = loss_shape, 
               aes(x = min_displacement, xend = max_displacement, 
                   y = percent_habitat_loss, yend = percent_habitat_loss,
                   group = scenario),
               color = "#3399FF", 
               size = 0.5, alpha = 0.8) +
  geom_point(data = loss_shape,
             aes(x = avg_displacement, y = percent_habitat_loss,
                 group = scenario),
             color = "#3399FF", 
             size = 2, alpha = 0.8) +
  
  # Location group 
  geom_segment(data = loss_loc, 
               aes(x = min_displacement, xend = max_displacement, 
                   y = percent_habitat_loss, yend = percent_habitat_loss,
                   group = location),
               color = "#109618", 
               size = 0.5, alpha = 0.8,
               position = position_jitter(height = 0.8, seed = 123)) +
  geom_point(data = loss_loc,
             aes(x = avg_displacement, y = percent_habitat_loss,
                 group = location),
             color = "#109618", 
             size = 2, alpha = 0.8,
             position = position_jitter(height = 0.8, seed = 123)) +
  
  # Add a manual legend
  annotate("point", x = 18, y = 84, color = "#FF9933", size = 2) +
  annotate("text", x = 20, y = 84, label = "Magnitude", hjust = 0) +
  annotate("point", x = 18, y = 77, color = "#3399FF", size = 2) +
  annotate("text", x = 20, y = 77, label = "Spatial var.", hjust = 0) +
  annotate("point", x = 18, y = 70, color = "#CC0000", size = 2) +
  annotate("text", x = 20, y = 70, label = "Size", hjust = 0) +
  annotate("point", x = 18, y = 63, color = "#109618", size = 2) +
  annotate("text", x = 20, y = 63, label = "Location", hjust = 0) +
  
  # Rest of plot specifications
  scale_y_continuous(limits=c(-1,102), breaks=(seq(0, 100, 20)), expand=c(0,0)) +
  scale_x_continuous(limits=c(0,33), breaks=(seq(0, 30, 10)), expand=c(0,0)) +
  labs(x = "Distance to analog (degrees)",
       y = "Species with complete\n habitat loss (%)") +
  theme_classic()






### Functions for 2d version

# Create the temperature gradient 
create_temperature_gradient <- function(grid_size, max_temp, min_temp) {
  # Create a SpatRaster with the specified dimensions
  r <- rast(ncol=grid_size, nrow=grid_size, xmin=0, xmax=grid_size, ymin=0, ymax=grid_size)
  
  # Create x and y coordinates for the entire grid
  xy <- crds(r)
  y <- xy[,2]
  
  # Calculate the middle y-coordinate (equator)
  middle_y <- grid_size/2
  
  # Calculate distance from equator for each cell (only considering y-axis)
  dist_from_equator <- abs(y - middle_y)
  
  # Maximum distance possible from equator to edge
  max_dist_from_equator <- grid_size/2
  
  # Normalize distance (0 at middle, 1 at extremes)
  normalized_dist <- dist_from_equator / max_dist_from_equator
  
  # Using the same quadratic function as temperature_gradient
  # that equals 1 at x=0 and 0 at x=1
  temp_factor <- 1 - normalized_dist^2
  
  # Calculate temperature (max_temp along equator, min_temp at top/bottom edges)
  temp <- min_temp + (max_temp - min_temp) * temp_factor
  
  # Assign values to the raster
  values(r) <- temp
  
  return(r)
}

# Generate random polygons of different sizes
generate_polygons <- function(num_polygons, grid_size, min_size, max_size = grid_size * 0.98) {
  # Create lists to store polygon data
  all_polygons <- list()
  ids <- numeric(num_polygons)
  
  # Center longitude for all polygons
  center_x <- grid_size / 2
  
  # Use the exact same distribution as in generate_segments "lognormal" option
  # beta distribution with alpha=0.3, beta=1
  alpha <- 0.3
  beta <- 1
  
  # Create positions array and apply the beta distribution
  positions <- seq(0, 1, length.out = num_polygons)
  size_ratios <- qbeta(positions, alpha, beta)
  
  # Sort to ensure a smooth distribution
  size_ratios <- sort(size_ratios)
  
  # Generate all polygons
  for (i in 1:num_polygons) {
    # Calculate target area as percentage of grid
    # Map size_ratio to area percentage (0-100%)
    area_pct <- size_ratios[i] * 100
    
    # Convert area percentage to dimensions
    # For very small areas, use circles
    if (area_pct < 15) {
      # Calculate radius from area
      area <- area_pct / 100 * grid_size * grid_size
      radius <- sqrt(area / pi)
      
      # Random y position within valid range
      pos_x <- center_x
      pos_y <- runif(1, radius, grid_size - radius)
      
      # Create circle-like polygon
      num_vertices <- max(6, round(8 * area_pct / 100 + 6))
      angles <- seq(0, 2*pi, length.out = num_vertices + 1)[1:num_vertices]
      
      # Generate coordinates
      x_coords <- pos_x + radius * cos(angles)
      y_coords <- pos_y + radius * sin(angles)
      
      # Close polygon
      x_coords <- c(x_coords, x_coords[1])
      y_coords <- c(y_coords, y_coords[1])
      
    } else {
      # For larger areas, use rectangles with fixed center longitude
      # Calculate height and width for target area
      
      # Height as a function of area percentage
      height_ratio <- min(0.98, 0.3 + (0.68 * area_pct / 100))
      height <- height_ratio * grid_size
      
      # Width calculated to achieve target area
      width <- (area_pct / 100 * grid_size * grid_size) / height
      
      # Center the rectangle along the center longitude
      x_min <- center_x - (width / 2)
      x_max <- center_x + (width / 2)
      
      # Position vertically with some variation for smaller polygons
      if (area_pct < 50) {
        max_offset <- (grid_size - height) / 2
        offset <- runif(1, -max_offset, max_offset) * (1 - area_pct/50)
        y_min <- (grid_size - height) / 2 + offset
      } else {
        # Center larger rectangles
        y_min <- (grid_size - height) / 2
      }
      
      y_max <- y_min + height
      
      # Create rectangle coordinates
      x_coords <- c(x_min, x_max, x_max, x_min, x_min)
      y_coords <- c(y_min, y_min, y_max, y_max, y_min)
    }
    
    # Ensure coordinates are within grid
    x_coords <- pmin(pmax(x_coords, 0), grid_size)
    y_coords <- pmin(pmax(y_coords, 0), grid_size)
    
    coords <- cbind(x_coords, y_coords)
    all_polygons[[i]] <- list(coords)
    ids[i] <- i
  }
  
  # Create SpatVector with all polygons
  polygons <- vect(all_polygons, type="polygons", crs="")
  polygons$id <- ids
  
  return(polygons)
} 

# Measure temperature for each polygon
measure_polygon_temperatures <- function(polygons, temp_raster) {
  # Extract values from raster for each polygon
  poly_temps <- extract(temp_raster, polygons)
  
  # Initialize results dataframe
  results <- data.frame(
    polygon_id = 1:length(polygons),
    area = expanse(polygons),
    mean_temp = NA,
    min_temp = NA,
    max_temp = NA,
    num_cells = NA
  )
  
  # Aggregate statistics by polygon ID
  poly_stats <- aggregate(poly_temps[, 2], by=list(ID=poly_temps$ID), 
                          FUN=function(x) c(mean=mean(x), min=min(x), max=max(x), count=length(x)))
  
  # Fill in results
  for (i in 1:nrow(poly_stats)) {
    id <- poly_stats[i, "ID"]
    stats <- poly_stats[i, "x"]
    
    results$mean_temp[id] <- stats[1]  # mean
    results$min_temp[id] <- stats[2]   # min
    results$max_temp[id] <- stats[3]   # max
    results$num_cells[id] <- stats[4]  # count
  }
  
  return(results)
}

# Measure temperature and exceeded envelopes for each polygon
measure_advanced_temperature <- function(polygons, current_raster, future_raster = NULL) {
  # Extract values from current raster for each polygon
  poly_temps_current <- extract(current_raster, polygons)
  
  # Initialize results dataframe
  results <- data.frame(
    polygon_id = 1:length(polygons),
    area = expanse(polygons),
    mean_temp = NA,
    min_temp = NA,
    max_temp = NA,
    num_cells = NA
  )
  
  # Aggregate statistics by polygon ID for current temperatures
  poly_stats <- aggregate(poly_temps_current[, 2], by=list(ID=poly_temps_current$ID), 
                          FUN=function(x) c(mean=mean(x), min=min(x), max=max(x), count=length(x)))
  
  # Fill in current temperature results
  for (i in 1:nrow(poly_stats)) {
    id <- poly_stats[i, "ID"]
    stats <- poly_stats[i, "x"]
    
    results$mean_temp[id] <- stats[1]  # mean
    results$min_temp[id] <- stats[2]   # min
    results$max_temp[id] <- stats[3]   # max
    results$num_cells[id] <- stats[4]  # count
  }
  
  # If future raster is provided, calculate additional metrics
  if (!is.null(future_raster)) {
    # Extract values from future raster for each polygon
    poly_temps_future <- extract(future_raster, polygons)
    
    # Add future temperature columns
    results$future_mean_temp <- NA
    results$future_min_temp <- NA
    results$future_max_temp <- NA
    results$temp_increase <- NA
    results$exceed_max_pct <- NA
    
    # Aggregate statistics for future temperatures
    future_poly_stats <- aggregate(poly_temps_future[, 2], by=list(ID=poly_temps_future$ID), 
                                   FUN=function(x) c(mean=mean(x), min=min(x), max=max(x)))
    
    # Fill in future temperature results
    for (i in 1:nrow(future_poly_stats)) {
      id <- future_poly_stats[i, "ID"]
      stats <- future_poly_stats[i, "x"]
      
      results$future_mean_temp[id] <- stats[1]  # mean
      results$future_min_temp[id] <- stats[2]   # min
      results$future_max_temp[id] <- stats[3]   # max
      
      # Calculate temperature increase
      results$temp_increase[id] <- results$future_mean_temp[id] - results$mean_temp[id]
    }
    
    # Calculate area exceeding maximum envelope for each polygon
    for (i in 1:length(polygons)) {
      # Check if the polygon has valid temperature data
      if (is.na(results$max_temp[i])) {
        results$exceed_max_pct[i] <- NA
        next
      }
      
      # Get current and future temperatures for this polygon
      poly_current <- poly_temps_current[poly_temps_current$ID == i, 2]
      poly_future <- poly_temps_future[poly_temps_future$ID == i, 2]
      
      # Ensure we have matching cells
      if (length(poly_current) != length(poly_future)) {
        warning(paste("Mismatch in cell counts for polygon", i))
        results$exceed_max_pct[i] <- NA
        next
      }
      
      # Get max temperature for this polygon from current climate
      threshold_temp <- results$max_temp[i]
      
      # Calculate percentage of cells that exceed the threshold
      if (length(poly_future) > 0) {
        pct_above_threshold <- sum(poly_future > threshold_temp) / length(poly_future) * 100
        results$exceed_max_pct[i] <- pct_above_threshold
      } else {
        results$exceed_max_pct[i] <- NA
      }
    }
  }
  
  return(results)
}

# Create sample segments representing ranges
create_latitude_segments <- function(
    latitude_sample,  # The latitude_sample data frame from your existing code
    segments_per_row = c(12, 6, 3, 1),  # Number of segments for each row
    y_spacing = 5,    # Vertical spacing between rows
    bottom_margin = -10, # Position of the bottom-most row
    gap_percentage = 15  # Percentage of total range to allocate to gaps
) {
  # Get latitude range from data
  latitude_min <- min(latitude_sample$y)
  latitude_max <- max(latitude_sample$y)
  total_range <- latitude_max - latitude_min
  
  # Split the data by scenario
  today_data <- latitude_sample %>% filter(scenario == "Today")
  future_data <- latitude_sample %>% filter(scenario == "Future")
  
  # Initialize output dataframe for segments
  all_segments <- data.frame()
  
  # Initialize dataframe for row statistics
  row_stats <- data.frame(
    row_id = integer(),
    avg_temp_increase = numeric(),
    pct_exceeding_envelope = numeric(),
    sd_pct_exceeding = numeric()
  )
  
  # Generate segments for each row
  for (row in 1:length(segments_per_row)) {
    num_segments <- segments_per_row[row]
    
    # Position for this row
    y_position <- bottom_margin - ((row - 1) * y_spacing)
    
    # Array to store percentages for each segment
    segment_percentages <- c()
    
    # Track statistics for this row
    row_segments <- data.frame()
    total_width <- 0
    weighted_temp_increase <- 0
    total_exceeding_width <- 0
    
    # Special handling for the last row (1 segment)
    if (num_segments == 1) {
      # For single segment, use full width without gaps
      segment_ranges <- list(c(latitude_min, latitude_max))
    } else {
      # Calculate how much space to allocate to segments vs. gaps
      total_gap_space <- (gap_percentage / 100) * total_range
      total_segment_space <- total_range - total_gap_space
      
      # Calculate individual segment width
      segment_width <- total_segment_space / num_segments
      
      # Calculate individual gap width
      gap_width <- total_gap_space / (num_segments - 1)
      
      # Create segments with specified widths and gaps
      segment_ranges <- list()
      current_start <- latitude_min
      
      for (i in 1:num_segments) {
        segment_end <- current_start + segment_width
        
        # Add this segment to our list
        segment_ranges[[i]] <- c(current_start, segment_end)
        
        # Move to next segment start (add gap)
        current_start <- segment_end + gap_width
      }
    }
    
    # Process each segment
    for (i in 1:length(segment_ranges)) {
      range_min <- segment_ranges[[i]][1]
      range_max <- segment_ranges[[i]][2]
      
      # For envelope determination, use dense sampling for accurate checking
      sample_points <- 100
      latitude_points <- seq(range_min, range_max, length.out = sample_points)
      
      # Get temperatures at these points
      today_temps <- approx(today_data$y, today_data$temperature, latitude_points)$y
      future_temps <- approx(future_data$y, future_data$temperature, latitude_points)$y
      
      # Calculate temperature increase for this segment
      temp_increase <- mean(future_temps - today_temps, na.rm = TRUE)
      
      # Find min/max of today's temperatures for this range
      min_envelope <- min(today_temps, na.rm = TRUE)
      max_envelope <- max(today_temps, na.rm = TRUE)
      
      # Check if future temperatures are within envelope
      within_envelope <- future_temps >= min_envelope & future_temps <= max_envelope
      
      # Calculate percentage exceeding envelope
      pct_exceeding <- mean(!within_envelope, na.rm = TRUE) * 100
      
      # Store this segment's percentage
      segment_percentages <- c(segment_percentages, pct_exceeding)
      
      # Use run-length encoding to find status changes
      rle_result <- rle(within_envelope)
      end_positions <- cumsum(rle_result$lengths)
      start_positions <- c(1, end_positions[-length(end_positions)] + 1)
      
      # Find which sample points correspond to status changes
      change_indices <- end_positions
      
      # Get the latitudes where status changes
      change_latitudes <- latitude_points[change_indices]
      
      # Add the start of range as first transition
      transition_points <- c(range_min, change_latitudes, range_max)
      
      # Create segments
      for (j in 1:(length(transition_points)-1)) {
        segment_width <- transition_points[j+1] - transition_points[j]
        
        # Fix: Properly determine the within_envelope status for this segment
        if (j <= length(rle_result$values)) {
          is_within <- rle_result$values[j]
        } else {
          # If we're past the end of the RLE values, use the last value
          is_within <- rle_result$values[length(rle_result$values)]
        }
        
        # Create segment from current transition to next
        segment_data <- data.frame(
          row_id = row,
          segment_id = paste(row, i, j, sep = "-"),
          y_start = transition_points[j],
          y_end = transition_points[j+1],
          y_position = y_position,
          width = segment_width,
          within_envelope = is_within,
          temp_increase = temp_increase
        )
        
        row_segments <- rbind(row_segments, segment_data)
        all_segments <- rbind(all_segments, segment_data)
        
        # Update statistics
        total_width <- total_width + segment_width
        weighted_temp_increase <- weighted_temp_increase + (segment_width * temp_increase)
        
        # Fix: Check for NA in is_within before using in the if condition
        if (!is.na(is_within) && !is_within) {
          total_exceeding_width <- total_exceeding_width + segment_width
        }
      }
    }
    
    # Calculate row statistics - avoid division by zero
    if (total_width > 0) {
      avg_temp_increase <- weighted_temp_increase / total_width
      pct_exceeding_envelope <- (total_exceeding_width / total_width) * 100
    } else {
      avg_temp_increase <- NA
      pct_exceeding_envelope <- NA
    }
    
    # Handle case where segment_percentages might be empty
    if (length(segment_percentages) > 0) {
      sd_pct_exceeding <- sd(segment_percentages, na.rm = TRUE)
    } else {
      sd_pct_exceeding <- NA
    }
    
    # Add row statistics
    row_stats <- rbind(row_stats, data.frame(
      row_id = row,
      avg_temp_increase = avg_temp_increase,
      pct_exceeding_envelope = pct_exceeding_envelope,
      sd_pct_exceeding = sd_pct_exceeding
    ))
  }
  
  # Return both segments and statistics
  return(list(
    segments = all_segments,
    stats = row_stats
  ))
}






# Run simulation 2d
# --------------

options(scipen=999)

# Parameters
grid_size <- 1000  # Grid dimensions (1000 x 1000)
num_polygons <- 1000  # Number of polygons to generate
min_polygon_size <- 2  # Minimum polygon size (radius for circular approximation)
max_polygon_size <- 1000  # Maximum polygon size

# Set random seed for reproducibility
set.seed(123)

# Create temperature gradient
temp_raster <- create_temperature_gradient(grid_size,max_temp=30, min_temp=0)

# Generate polygons
polygons <- generate_polygons(num_polygons, grid_size, min_polygon_size, max_polygon_size)

hist(100 * expanse(polygons) / (grid_size * grid_size), breaks = seq(0, 100, by = 5), 
     main = "Polygon Areas", xlab = "Area (% of grid)", col = "lightgray")

# Create future temperature gradient
future_temp_raster <- create_temperature_gradient(grid_size, max_temp=33, min_temp=3)

# Calculate temperature increase and envelope exceedance
results_2d <- measure_advanced_temperature(polygons, temp_raster, future_temp_raster)



## plot latitude/longitude and temperature profile in 2d

temp_df <- as.data.frame(temp_raster, xy=TRUE)
names(temp_df) <- c("x", "y", "temperature")
temp_df$scenario <- "Today"

future_temp_df <- as.data.frame(future_temp_raster, xy=TRUE)
names(future_temp_df) <- c("x", "y", "temperature")
future_temp_df$scenario <- "Future"
combined_temp_df <- rbind(temp_df, future_temp_df)

latitude_sample <- combined_temp_df[combined_temp_df$x == 0.5, ]
plot1_2d <- ggplot(latitude_sample, aes(x = y, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Temperature (°C)",
       color = "Scenario")+
  scale_y_continuous(limits = c(0,36),breaks=c(seq(0,35,5)), expand = c(0,0))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


segment_result <- create_latitude_segments(latitude_sample)
segments_plot <- segment_result$segments
row_stats <- segment_result$stats

y_position_map <- c(
  "-10" = -4, 
  "-15" = -5,    
  "-20" = -6,  
  "-25" = -7    
)

# Apply the mapping to create a new y_position column
segments_plot$y_position_new <- y_position_map[as.character(segments_plot$y_position)]

# Format the statistics for display
stats_text <- row_stats %>%
  mutate(
    label = sprintf("%.1f°C, %.1f%% ± %.1f%%",
                    avg_temp_increase, 
                    pct_exceeding_envelope,
                    sd_pct_exceeding))


plot_fig1 <- ggplot(combined_temp, aes(x = latitude, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Temperature (°C)",
       color = "Scenario")+
  scale_x_continuous(breaks = seq(-50, 50, by = 50)) 
#theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())


plot_fig1 + 
  geom_segment(data = segments_plot,
               aes(x = y_start, xend = y_end,
                   y = y_position_new, yend = y_position_new,
                   color = within_envelope),
               linewidth = 1.5) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red", 
                                "TRUE" = "#3d8c40", "FALSE" = "orange"),
                     name = "Temperature/Status",
                     breaks = c("Today", "Future", TRUE, FALSE),
                     labels = c("Today", "Future", "Within envelope", "Exceeding envelope")) +
  # Add statistics text annotations
  # geom_text(data = stats_text,
  #           aes(x = max(segments_plot$y_end) + 20, # Position to the right of the segments
  #               y = unique(segments_plot$y_position)[row_id], # Position at each row's y value
  #               label = label),
  #           hjust = 0, # Left-align text
  #           size = 3.5) +
  # Set y-axis to start at 0
  scale_y_continuous(
    limits = c(min(segments_plot$y_position_new) - 0.5, NA),
    breaks = seq(0, 40, by = 10),  # Only show positive breaks
    labels = function(x) ifelse(x >= 0, x, "")  # Hide negative labels
  ) +
  # Add a horizontal line at y=0 to serve as x-axis
  geom_hline(yintercept = 0) +
  # Remove all x-axis elements from the theme
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),  # Remove the default x title
    axis.line.x = element_blank(),
    # Fix legend to show only lines
    legend.key.size = unit(1.5, "lines"),
    legend.key.height = unit(1.5, "lines")
  ) +
  # Add a single custom "Latitude" label at y=0
  annotate("text", x = 500, y = 0, 
           label = "Latitude", vjust = 2.5, hjust = 0.5) +
  # Remove the original x-axis title
  xlab(NULL) +
  # Make sure plot extends to show the text
  coord_cartesian(xlim = c(NA, max(segments_plot$y_end) + 200)) +
  # Override the line segments in the legend to show only lines
  guides(color = guide_legend(override.aes = list(shape = NA)))






#linear version

# Create x, y values
x <- rep(0.5, 1000 * 2)
y <- rep(seq(0.5, 999.5, by = 1), 2)

# Create temperature values (linear relationship)
# Today: 0 to 30
# Future: 3 to 33
temp_today_lin <- seq(0, 30, length.out = 1000)
temp_future_lin <- seq(3, 33, length.out = 1000)
temperature <- c(temp_today_lin, temp_future_lin)

# Create scenario values
scenario <- rep(c("Today", "Future"), each = 1000)

# Combine into a dataframe
latitude_sample <- data.frame(
  x = x,
  y = y,
  temperature = temperature,
  scenario = scenario)


segment_result <- create_latitude_segments(latitude_sample)
segments_plot <- segment_result$segments
row_stats <- segment_result$stats

y_position_map <- c(
  "-10" = -4, 
  "-15" = -5,    
  "-20" = -6,  
  "-25" = -7    
)

# Apply the mapping to create a new y_position column
segments_plot$y_position_new <- y_position_map[as.character(segments_plot$y_position)]

# Format the statistics for display
stats_text <- row_stats %>%
  mutate(
    label = sprintf("%.1f°C, %.1f%% ± %.1f%%",
                    avg_temp_increase, 
                    pct_exceeding_envelope,
                    sd_pct_exceeding))


plot_fig1 <- ggplot(latitude_sample, aes(x = y, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Temperature (°C)",
       color = "Scenario")+
  scale_x_continuous(breaks = seq(-50, 50, by = 50)) 
#theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())


plot_fig1 + 
  geom_segment(data = segments_plot,
               aes(x = y_start, xend = y_end,
                   y = y_position_new, yend = y_position_new,
                   color = within_envelope),
               linewidth = 1.5) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red", 
                                "TRUE" = "#3d8c40", "FALSE" = "orange"),
                     name = "Temperature/Status",
                     breaks = c("Today", "Future", TRUE, FALSE),
                     labels = c("Today", "Future", "Within envelope", "Exceeding envelope")) +
  # Add statistics text annotations
  # geom_text(data = stats_text,
  #           aes(x = max(segments_plot$y_end) + 20, # Position to the right of the segments
  #               y = unique(segments_plot$y_position)[row_id], # Position at each row's y value
  #               label = label),
  #           hjust = 0, # Left-align text
  #           size = 3.5) +
  # Set y-axis to start at 0
  scale_y_continuous(
    limits = c(min(segments_plot$y_position_new) - 0.5, NA),
    breaks = seq(0, 40, by = 10),  # Only show positive breaks
    labels = function(x) ifelse(x >= 0, x, "")  # Hide negative labels
  ) +
  # Add a horizontal line at y=0 to serve as x-axis
  geom_hline(yintercept = 0) +
  # Remove all x-axis elements from the theme
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),  # Remove the default x title
    axis.line.x = element_blank(),
    # Fix legend to show only lines
    legend.key.size = unit(1.5, "lines"),
    legend.key.height = unit(1.5, "lines")
  ) +
  # Add a single custom "Latitude" label at y=0
  annotate("text", x = 500, y = 0, 
           label = "Latitude", vjust = 2.5, hjust = 0.5) +
  # Remove the original x-axis title
  xlab(NULL) +
  # Make sure plot extends to show the text
  coord_cartesian(xlim = c(NA, max(segments_plot$y_end) + 200)) +
  # Override the line segments in the legend to show only lines
  guides(color = guide_legend(override.aes = list(shape = NA)))






lin_lat_today <- latitude_sample[latitude_sample$scenario == "Today", ]
names(lin_lat_today)[names(lin_lat_today) == "y"] <- "latitude"
lin_lat_future <- latitude_sample[latitude_sample$scenario == "Future", ]
names(lin_lat_future)[names(lin_lat_future) == "y"] <- "latitude"
linear_displacement <- calculate_latitudinal_displacement(lin_lat_today, lin_lat_future)








# Plot displacement vs. latitude
ggplot() +
  geom_line(data=linear_displacement, aes(x = display_latitude, y = display_displacement),size = 1,linetype="dashed") +
  geom_line(data=displacement, aes(x = display_latitude, y = display_displacement),size = 1) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Displacement \n(degrees latitude)")+
  scale_x_continuous(limits = c(-50, 50),expand=c(0,0))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))













plot2_2d <- ggplot(results_2d, aes(x = area/ (grid_size * grid_size)*100, y = temp_increase)) +
  geom_pointdensity(position = position_jitter(width = 0.01, height = 0.05),size=0.01)+
  scale_color_gradient(low='#e9f6fd',high="#2171b5") + 
  theme_classic() +
  labs(x = "Area (% \nof domain)", 
       y = "Temperature \nexposure (°C)")+
  scale_x_continuous(limits = c(0, 100), 
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),)+
  theme(legend.position="none")

results_2d$pct_domain <- results_2d$area / (grid_size*grid_size)

results_2d <- results_2d %>%
  arrange(pct_domain) %>%
  mutate(cum_species_pct = (row_number() / n()) * 100)

# plot range size vs habitat loss
plot3_2d <- ggplot(results_2d, aes(x = pct_domain*100, y = exceed_max_pct)) +
  geom_pointdensity(size=0.8)+
  scale_color_gradient(low='#e9f6fd',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  geom_line(aes(x = pct_domain*100,y = cum_species_pct * max(exceed_max_pct) / 100),
            color = "blue", linewidth = 1) +
  scale_y_continuous(
    name = "Area exceeding\n envelope (%)",
    sec.axis = sec_axis(~. * 100 / max(results_sim$exceed_max_pct), 
                        name = "Cumulative species (%)",breaks = seq(0, 100, 20)),
    limits = c(0, 100.5),
    breaks = seq(0, 100, 20),
    expand=c(0,0)) +
  theme_classic() +
  labs(x = "Proportion \nof domain (%)", 
       y = "Area exceeding envelope (%)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  theme(legend.position="none",
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.line.y.right = element_line(color="blue"))

#plot range size distribution
plot4_2d <- ggplot(results_2d,aes(x=area / (grid_size * grid_size)))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Area (% \nof domain)", 
       y = "Frequency \nof species")+
  scale_x_continuous(limits = c(-0.05, 1))


#plot temperature raster with example polygons
areas <- expanse(polygons)
polygons$area <- areas
sorted_polys <- polygons[order(polygons$area),]
selected_indices <- round(seq(1, length(sorted_polys), length.out = 5))
selected_polys <- sorted_polys[selected_indices,]
selected_polys <- selected_polys[1:3,]
selected_polys_sf <- st_as_sf(selected_polys)
selected_polys_sf$poly_id <- 1:nrow(selected_polys_sf)
r <- raster(temp_raster)
rdf <- as.data.frame(r, xy = TRUE)
colnames(rdf) <- c("x", "y", "temp")

temprast_poly <- ggplot() +
  geom_raster(data = rdf, aes(x = x, y = y, fill = temp)) +
  scale_fill_gradientn(
    colors = hcl.colors(100, "RdYlBu", rev = TRUE),
    name = "Temperature",
    na.value = NA) +
  geom_sf(data = selected_polys_sf, fill = NA, 
          color = "black", linewidth = 1) +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)) 





sim_rastpoly <- plot_grid(temprast_poly,labels = c("A"), ncol = 2)

title_2d <- ggdraw()+draw_label("2-D",fontface="bold",x=0.1,y=0.5)
sim_2d <- plot_grid(plot1_2d, plot2_2d,plot3_2d,plot4_2d, labels = c( "B","C","D","E"), ncol = 4)
sim_2dt <- plot_grid(title_2d,sim_2d,ncol=1,rel_heights=c(0.1,1))

title_1d <- ggdraw()+draw_label("1-D",fontface="bold",x=0.1,y=0.5)
sim_1d <- plot_grid(plot1, plot2,plot3_nored,plot4, labels = c( "F","G","H","I"), ncol = 4)
sim_1dt <- plot_grid(title_1d,sim_1d,ncol=1,rel_heights=c(0.1,1))


fig_s_2d <- plot_grid(sim_rastpoly, sim_2dt, sim_1dt, ncol = 1, rel_heights = c(1.2, 1, 1))


ggsave(paste0(figures_folder, "fig_S_2d.png"), fig_s_2d, 
       width = 11, height = 8, dpi = 350)























#  SSP temperature data scenario

# Create a smooth temperature function based on the SSP data
create_temp_gradient_from_data <- function(lat_temp_data, latitude_range, span = 0.5) {
  # Create a sequence of latitudes covering the full range
  lat_seq <- seq(latitude_range[1], latitude_range[2], length.out = 1000)
  
  # Calculate the midpoint in terms of the sequence index
  midpoint_index <- ceiling(length(lat_seq) / 2)
  midpoint_lat <- lat_seq[midpoint_index]
  
  # Normalize the actual latitude data to the range we're using
  min_actual_lat <- min(lat_temp_data$latitude)
  max_actual_lat <- max(lat_temp_data$latitude)
  actual_lat_range <- max_actual_lat - min_actual_lat
  
  # Sort the data by latitude to ensure proper ordering
  lat_temp_data <- lat_temp_data[order(lat_temp_data$latitude),]
  
  # Get only the first half of the data (up to median latitude)
  median_lat <- median(lat_temp_data$latitude)
  first_half <- lat_temp_data[lat_temp_data$latitude <= median_lat,]
  
  # If we have too few points, include some from just past the median
  if(nrow(first_half) < 5) {
    n_add <- min(5 - nrow(first_half), nrow(lat_temp_data) - nrow(first_half))
    extra_points <- lat_temp_data[order(lat_temp_data$latitude)][nrow(first_half) + (1:n_add),]
    first_half <- rbind(first_half, extra_points)
  }
  
  # Fit a loess curve to the first half
  # Using a higher span for smoother curve
  loess_fit <- loess(temperature ~ latitude, data = first_half, span = span)
  
  # Create temperature array
  temps <- numeric(length(lat_seq))
  
  # Calculate temperatures for the first half (including midpoint)
  for (i in 1:midpoint_index) {
    # Convert simulation latitude to real-world latitude for lookup
    norm_pos <- (lat_seq[i] - latitude_range[1]) / (latitude_range[2] - latitude_range[1])
    real_lat <- min_actual_lat + norm_pos * actual_lat_range
    
    # Ensure we stay within the bounds of the real data
    real_lat <- max(min_actual_lat, min(max_actual_lat, real_lat))
    
    # Get temperature at this latitude using the loess fit
    temps[i] <- predict(loess_fit, real_lat)
  }
  
  # Check if there's a dip near the middle by looking at the slope
  # Calculate the slope for the last few points in the first half
  slopes <- diff(temps[(midpoint_index-10):midpoint_index]) / 
    diff(lat_seq[(midpoint_index-10):midpoint_index])
  
  # If the slope turns negative before the midpoint, adjust the curve
  if(any(slopes < 0)) {
    # Find the index where the slope starts to decrease
    max_temp_idx <- which.max(temps[1:midpoint_index])
    
    # Set all temperatures between max_temp_idx and midpoint_index to the maximum
    if(max_temp_idx < midpoint_index) {
      temps[max_temp_idx:midpoint_index] <- temps[max_temp_idx]
    }
  }
  
  # Mirror the first half for the second half
  for (i in (midpoint_index+1):length(lat_seq)) {
    # Find the corresponding index in the first half
    mirror_index <- 2 * midpoint_index - i
    # Make sure we don't go below the first index
    mirror_index <- max(1, mirror_index)
    temps[i] <- temps[mirror_index]
  }
  
  # Create a temperature gradient object in the format expected
  temp_gradient <- list(
    latitude = lat_seq,
    temperature = temps
  )
  
  return(temp_gradient)
} 

hist_today_data <- combined_lat585[combined_lat585$scenario == "Today",]
hist_today_data <- hist_today_data[order(today_data$latitude),]
hist_today_temp <- create_temp_gradient_from_data(hist_today_data, c(0, 1000))
hist_today_df <- data.frame(latitude = hist_today_temp$latitude,temperature = hist_today_temp$temperature,scenario = "Today")

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


sim_clim585 <- rbind(hist_today_df, ssp585_future_df)
sim_clim370 <- rbind(hist_today_df, ssp370_future_df)
sim_clim245 <- rbind(hist_today_df, ssp245_future_df)
sim_clim126 <- rbind(hist_today_df, ssp126_future_df)


ggplot(sim_clim585, aes(x = latitude, y = temperature, color = scenario)) +
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
ssp585_results <- analyze_temp_exceed(segments, hist_today_temp, ssp585_future_temp)
ssp585_results$scenario <- "SSP585"
ssp585_results$min_change <- NA  
ssp585_results$max_change <- NA  
ssp585_results$pct_domain <- results$size / 1000

ssp370_results <- analyze_temp_exceed(segments, hist_today_temp, ssp370_future_temp)
ssp370_results$scenario <- "SSP370"
ssp370_results$min_change <- NA  
ssp370_results$max_change <- NA  
ssp370_results$pct_domain <- results$size / 1000

ssp245_results <- analyze_temp_exceed(segments, hist_today_temp, ssp245_future_temp)
ssp245_results$scenario <- "SSP245"
ssp245_results$min_change <- NA  
ssp245_results$max_change <- NA  
ssp245_results$pct_domain <- results$size / 1000

ssp126_results <- analyze_temp_exceed(segments, hist_today_temp, ssp126_future_temp)
ssp126_results$scenario <- "SSP126"
ssp126_results$min_change <- NA  
ssp126_results$max_change <- NA  
ssp126_results$pct_domain <- results$size / 1000

# Combine with other results

ssp_results <- rbind(ssp585_results, ssp370_results,ssp245_results,ssp126_results)

#add ssp scenarios to tempincr
ssp585_increase <- data.frame(
  latitude = hist_today_temp$latitude,
  temp_increase = ssp585_future_temp$temperature - hist_today_temp$temperature,
  scenario = "SSP585")

ssp370_increase <- data.frame(
  latitude = hist_today_temp$latitude,
  temp_increase = ssp370_future_temp$temperature - hist_today_temp$temperature,
  scenario = "SSP370")

ssp245_increase <- data.frame(
  latitude = hist_today_temp$latitude,
  temp_increase = ssp245_future_temp$temperature - hist_today_temp$temperature,
  scenario = "SSP245")

ssp126_increase <- data.frame(
  latitude = hist_today_temp$latitude,
  temp_increase = ssp126_future_temp$temperature - hist_today_temp$temperature,
  scenario = "SSP126")

# Combine all temperature increase data
ssp_temp_increases <- rbind(ssp585_increase, ssp370_increase, ssp245_increase, ssp126_increase)

ssp_temp_increases$latitude_degrees <- -50 + 100 * (ssp_temp_increases$latitude - min(ssp_temp_increases$latitude)) / 
  (max(ssp_temp_increases$latitude) - min(ssp_temp_increases$latitude))


displacement_585 <- calculate_latitudinal_displacement(hist_today_df, ssp585_future_df)
displacement_370 <- calculate_latitudinal_displacement(hist_today_df, ssp370_future_df)
displacement_245 <- calculate_latitudinal_displacement(hist_today_df, ssp245_future_df)
displacement_126 <- calculate_latitudinal_displacement(hist_today_df, ssp126_future_df)

analyzed_segments_585 <- find_critical_range_size(displacement_585, segments)
analyzed_segments_370 <- find_critical_range_size(displacement_370, segments)
analyzed_segments_245 <- find_critical_range_size(displacement_245, segments)
analyzed_segments_126 <- find_critical_range_size(displacement_126, segments)

# Add  identifier to both datasets
displacement_585$scenario <- "ssp585"
analyzed_segments_585$scenario <- "ssp585"
displacement_370$scenario <- "ssp370"
analyzed_segments_370$scenario <- "ssp370"
displacement_245$scenario <- "ssp245"
analyzed_segments_245$scenario <- "ssp245"
displacement_126$scenario <- "ssp126"
analyzed_segments_126$scenario <- "ssp126"

# Combine with previous results
displacement_by_ssp <- rbind(displacement_585,displacement_370,displacement_245,displacement_126)
analyzed_segments_by_ssp <- rbind(analyzed_segments_585,analyzed_segments_370,analyzed_segments_245,analyzed_segments_126)


# Set a custom color palette
scenario_colors <- c("SSP585" = "black","SSP370" = "darkgrey","SSP245" = "grey","SSP126"="lightgrey")

# Plot for habitat loss vs range size
comparison_plot <- ggplot(ssp_results, aes(x = pct_domain*100, y = exceed_max_pct, color = scenario)) +
  #stat_summary_bin(fun = "mean", geom = "line",linewidth = 1.2,bins = 10) +
  geom_smooth(method = "loess", span = 0.3,se=F,linewidth=1.2) +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Habitat loss\n(% above today's max temp)",
       color = "Warming\nScenario") +
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0),
                     limits=c(0,100))+
  scale_color_manual(values = scenario_colors) +
  theme(legend.position = "right")+
  theme(legend.position = "none")

# Plot for mean temperature increase vs range size
plot2 <-  ggplot(ssp_temp_increases, aes(x = latitude_degrees, y = temp_increase, color = scenario)) +
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
  scale_x_continuous(breaks = seq(-50, 50, by = 50),expand=c(0,0))+
  scale_y_continuous(limits=c(0,4),expand=c(0,0))+
  theme(legend.position = "none")


legend_plot <- ggplot(ssp_temp_increases, aes(x = latitude_degrees, y = temp_increase, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = scenario_colors,
                     name = "Scenario") +
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

# Plot Latitudinal displacement plot
lat_disp_ssp <- ggplot(analyzed_segments_by_ssp, aes(x = display_center_lat, y = display_size)) +
  # Regular points
  geom_point(color = "black", size = 0.5, alpha = 0.1) +
  # No analog points 
  geom_point(data = subset(analyzed_segments_by_ssp, no_analog), 
             aes(x = display_center_lat, y = display_size),
             color = "black", size = 0.5, alpha = 0.1) +
  # required latitudinal shift
  geom_line(data = displacement_by_ssp, 
            aes(x = display_latitude, y = abs(display_displacement), 
                color = factor(scenario), group = factor(scenario)), 
            size = 1.2) +
  scale_y_continuous(limits=c(-2,100),expand=c(0,0))+
  scale_x_continuous(limits=c(-50,50),expand=c(0,0))+
  scale_color_manual(values = c("lightgrey", "grey","darkgrey","black")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  theme(legend.position = "none")

complete_loss_ssp <- analyzed_segments_by_ssp %>%
  group_by(scenario) %>%
  summarize(
    total_count = n(),
    habitat_loss_count = sum(complete_habitat_loss),
    percentage = (habitat_loss_count / total_count) * 100)

complete_loss_ssp_plot <- ggplot(complete_loss_ssp, aes(x = scenario, y = percentage, fill = scenario)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("lightgrey", "grey","darkgrey","black")) +  
  labs(x = "Scenario",
       y = "Species with \ncomplete habitat loss (%)",
       fill = "Scenario") +
  scale_y_continuous(limits=c(0,100), expand=c(0,0))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

sim_ssp <- plot_grid(plot2, comparison_plot, lat_disp_ssp,complete_loss_ssp_plot,
                     legend,labels = c("A", "B","C","D",""), ncol = 5,rel_widths = c(1, 1, 1,1,0.4),align="hv",axis = "tblr")

sim_ssp585 <- create_plot_gridBC(hist_today_df%>% select(-scenario), ssp585_future_df%>% select(-scenario),
                                 ssp_results[ssp_results$scenario == "SSP585",],
                                 labels=c("E","F","G","H"), title="SSP585")
sim_ssp370 <- create_plot_gridBC(hist_today_df%>% select(-scenario), ssp370_future_df%>% select(-scenario),
                                 ssp_results[ssp_results$scenario =="SSP370",],
                                 labels=c("I","J","K","L"), title="SSP370")
sim_ssp245 <- create_plot_gridBC(hist_today_df%>% select(-scenario), ssp245_future_df%>% select(-scenario),
                                 ssp_results[ssp_results$scenario =="SSP245",],
                                 labels=c("M","N","O","P"), title="SSP245")
sim_ssp126 <- create_plot_gridBC(hist_today_df%>% select(-scenario), ssp126_future_df%>% select(-scenario),
                                 ssp_results[ssp_results$scenario =="SSP126",],
                                 labels=c("Q","R","S","T"), title="SSP126")

fig_s_ssp <- plot_grid(
  sim_ssp, sim_ssp585, sim_ssp370, sim_ssp245, sim_ssp126, 
  ncol = 1, rel_heights = c(1.3, 1, 1, 1, 1))


ggsave(paste0(figures_folder, "fig_S_ssp.png"), fig_s_ssp, 
       width = 12, height = 13, dpi = 350)







