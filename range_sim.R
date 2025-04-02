# range size simulation

library(terra)
library(ggplot2)
library(dplyr)
library(ggpointdensity)





# Create a temperature gradient across latitude
temperature_gradient <- function(max_temp = 30, min_temp = 0, latitude_range = c(0, 1000)) {
  # Create a sequence of latitudinal points
  min_lat <- latitude_range[1]
  max_lat <- latitude_range[2]
  latitude <- seq(min_lat, max_lat, by=0.5)
  
  # Middle latitude (equator)
  middle_lat <- (max_lat + min_lat)/2
  
  # Normalize distance from middle (0 at middle, 1 at extremes)
  dist_from_middle <- abs(latitude - middle_lat)
  max_dist <- max(dist_from_middle)  # This would be half the latitude range
  normalized_dist <- dist_from_middle / max_dist
  
  # The temperature should be max_temp at the middle and min_temp at the extremes
  # Using a simple quadratic function that equals 1 at x=0 and 0 at x=1
  temp_factor <- 1 - normalized_dist^2
  
  # Calculate temperature
  temp <- min_temp + (max_temp - min_temp) * temp_factor
  
  # Return data frame with latitude and temperature
  return(data.frame(latitude = latitude, temperature = temp))
}

# Generate random segments of different sizes along latitude
#distribution options: uniform normal left_skewed right_skewed lognormal inverse_lognormal
generate_segments <- function(num_segments, min_size, max_size, latitude_range = c(0, 1000), 
                              distribution = "uniform") {
  # Get the range of latitude
  min_lat <- latitude_range[1]
  max_lat <- latitude_range[2]
  lat_span <- max_lat - min_lat
  
  # Check if max_size is reasonable for the latitude range
  if (max_size > lat_span) {
    warning(paste("max_size (", max_size, ") is larger than the latitude span (", 
                  lat_span, "). Adjusting to 1/4 of the span."))
    max_size <- lat_span / 4
  }
  
  # Create data frame to store segment data
  segments <- data.frame(
    segment_id = 1:num_segments,
    start = numeric(num_segments),
    end = numeric(num_segments)
  )
  
  for (i in 1:num_segments) {
    valid_segment <- FALSE
    attempt <- 0
    max_attempts <- 100
    
    while (!valid_segment && attempt < max_attempts) {
      attempt <- attempt + 1
      
      # Generate random size within limits based on distribution
      size <- if (distribution == "uniform") {
        # Even distribution (uniform)
        round(runif(1, min_size, max_size))
      } else if (distribution == "normal") {
        # Normal distribution centered between min and max
        mean_size <- (min_size + max_size) / 2
        sd_size <- (max_size - min_size) / 6  # ~99.7% of values within range
        # Ensure the size stays within bounds
        round(max(min_size, min(max_size, rnorm(1, mean_size, sd_size))))
      } else if (distribution == "left_skewed") {
        # Left skewed (more larger sizes)
        # Using beta distribution with shape parameters favoring larger values
        round(min_size + (max_size - min_size) * rbeta(1, 5, 2))
      } else if (distribution == "right_skewed") {
        # Right skewed (more smaller sizes)
        # Using beta distribution with shape parameters favoring smaller values
        round(min_size + (max_size - min_size) * rbeta(1, 2, 5))
      } else if (distribution == "lognormal") {
        # Very right-skewed 
        # Using lognormal with sigma=2
        #meanlog <- log(min_size) # Set log-scale mean to create desired distribution
        #sigma <- 2
        # Generate from lognormal and cap at max_size
        #round(min(max_size, exp(rnorm(1, meanlog, sigma))))
        
        # Create a hollow curve distribution 
        alpha <- 0.3  # Very small alpha creates strong right skew (many small values)
        beta <- 1   # Higher beta concentrates values toward small end
        # Scale to desired range
        size <- min_size + (max_size - min_size) * rbeta(1, alpha, beta)
        round(size)
      } else if (distribution == "inverse_lognormal") {
        # Inverse of lognormal - many large sizes, few small ones
        # First generate lognormal
        #meanlog <- log(min_size)
        #sigma <- 2
        #log_val <- exp(rnorm(1, meanlog, sigma))
        # Then invert the distribution within our range
        #round(max(min_size, max_size - min(max_size - min_size, log_val)))
        
        # Inverse distribution - more large segments
        alpha <- 3.0 
        beta <- 0.3  # Very small beta creates strong left skew (many large values)
        # Scale to desired range
        size <- min_size + (max_size - min_size) * rbeta(1, alpha, beta)
        round(size)
      } else {
        # Default to uniform if invalid option is provided
        round(runif(1, min_size, max_size))
      }
      
      # Generate random center position
      center <- round(runif(1, min_lat + size/2, max_lat - size/2))
      
      # Calculate start and end points
      start <- center - size/2
      end <- center + size/2
      
      # Ensure segment is within bounds
      if (start >= min_lat && end <= max_lat) {
        valid_segment <- TRUE
        segments$start[i] <- start
        segments$end[i] <- end
      }
    }
    
    # If we couldn't find a valid segment after max attempts, place one at the center
    if (!valid_segment) {
      center <- (min_lat + max_lat) / 2
      size <- min(max_size, lat_span / 10)
      segments$start[i] <- center - size/2
      segments$end[i] <- center + size/2
      warning(paste("Couldn't find valid placement for segment", i, "after", max_attempts, "attempts. Placing at center."))
    }
  }
  
  return(segments)
}

# Generate random segments with smooth size distributions along latitude
#distribution options: uniform normal left_skewed right_skewed lognormal inverse_lognormal
generate_segments <- function(num_segments, min_size, max_size, latitude_range = c(0, 1000), 
                              distribution = "uniform") {
  # Get the range of latitude
  min_lat <- latitude_range[1]
  max_lat <- latitude_range[2]
  lat_span <- max_lat - min_lat
  
  # Check if max_size is reasonable for the latitude range
  if (max_size > lat_span) {
    warning(paste("max_size (", max_size, ") is larger than the latitude span (", 
                  lat_span, "). Adjusting to 1/4 of the span."))
    max_size <- lat_span / 4
  }
  
  # Create data frame to store segment data
  segments <- data.frame(
    segment_id = 1:num_segments,
    start = numeric(num_segments),
    end = numeric(num_segments)
  )
  
  # First, generate all segment sizes at once using a smooth distribution
  positions <- seq(0, 1, length.out = num_segments)
  
  # Generate sizes based on distribution type
  sizes <- if (distribution == "uniform") {
    # Evenly distributed sizes
    seq(min_size, max_size, length.out = num_segments)
  } else if (distribution == "normal") {
    # Normal distribution
    mean_size <- (min_size + max_size) / 2
    sd_size <- (max_size - min_size) / 6
    sapply(positions, function(p) {
      round(max(min_size, min(max_size, mean_size + sd_size * qnorm(p))))
    })
  } else if (distribution == "left_skewed") {
    # Left skewed (more larger sizes)
    sapply(positions, function(p) {
      round(min_size + (max_size - min_size) * qbeta(p, 5, 2))
    })
  } else if (distribution == "right_skewed") {
    # Right skewed (more smaller sizes)
    sapply(positions, function(p) {
      round(min_size + (max_size - min_size) * qbeta(p, 2, 5))
    })
  } else if (distribution == "lognormal") {
    # Very right-skewed distribution
    alpha <- 0.3
    beta <- 1
    sapply(positions, function(p) {
      round(min_size + (max_size - min_size) * qbeta(p, alpha, beta))
    })
  } else if (distribution == "inverse_lognormal") {
    # More large segments
    alpha <- 3.0
    beta <- 0.3
    sapply(positions, function(p) {
      round(min_size + (max_size - min_size) * qbeta(p, alpha, beta))
    })
  } else {
    # Default to uniform
    seq(min_size, max_size, length.out = num_segments)
  }
  
  # Now randomly place segments with these predefined sizes
  for (i in 1:num_segments) {
    valid_segment <- FALSE
    attempt <- 0
    max_attempts <- 100
    size <- sizes[i]
    
    while (!valid_segment && attempt < max_attempts) {
      attempt <- attempt + 1
      
      # Generate random center position
      center <- round(runif(1, min_lat + size/2, max_lat - size/2))
      
      # Calculate start and end points
      start <- center - size/2
      end <- center + size/2
      
      # Ensure segment is within bounds
      if (start >= min_lat && end <= max_lat) {
        valid_segment <- TRUE
        segments$start[i] <- start
        segments$end[i] <- end
      }
    }
    
    # If we couldn't find a valid segment after max attempts, place one at the center
    if (!valid_segment) {
      center <- (min_lat + max_lat) / 2
      size <- min(size, lat_span / 10)
      segments$start[i] <- center - size/2
      segments$end[i] <- center + size/2
      warning(paste("Couldn't find valid placement for segment", i, "after", max_attempts, "attempts. Placing at center."))
    }
  }
  
  return(segments)
}

# Measure temperature and exceeded envelopes for each segment
analyze_temp_exceed <- function(segments, current_gradient, future_gradient = NULL) {
  # Create results data frame
  results <- data.frame(
    segment_id = segments$segment_id,
    size = segments$end - segments$start,
    mean_temp = numeric(nrow(segments)),
    min_temp = numeric(nrow(segments)),
    max_temp = numeric(nrow(segments)))
  
  for (i in 1:nrow(segments)) {
    # Get start and end points
    start <- segments$start[i]
    end <- segments$end[i]
    
    # Find temperature points within this segment
    segment_indices <- which(current_gradient$latitude >= start & 
                               current_gradient$latitude <= end)
    
    # Get temperatures for this segment
    segment_temps <- current_gradient$temperature[segment_indices]
    
    # Calculate statistics
    results$mean_temp[i] <- mean(segment_temps)
    results$min_temp[i] <- min(segment_temps)
    results$max_temp[i] <- max(segment_temps)
  }
  
  # If future gradient is provided, calculate additional metrics
  if (!is.null(future_gradient)) {
    # Add future temperature columns
    results$future_mean_temp <- numeric(nrow(segments))
    results$future_min_temp <- numeric(nrow(segments))
    results$future_max_temp <- numeric(nrow(segments))
    results$temp_increase <- numeric(nrow(segments))
    results$exceed_max_pct <- numeric(nrow(segments))
    results$exceed_min_pct <- numeric(nrow(segments))
    results$exceed_envelope_pct <- numeric(nrow(segments))
    
    for (i in 1:nrow(segments)) {
      # Get start and end points
      start <- segments$start[i]
      end <- segments$end[i]
      
      # Find points within this segment
      segment_indices <- which(future_gradient$latitude >= start & 
                                 future_gradient$latitude <= end)
      
      # Get future temperatures for this segment
      future_temps <- future_gradient$temperature[segment_indices]
      
      # Calculate statistics
      results$future_mean_temp[i] <- mean(future_temps)
      results$future_min_temp[i] <- min(future_temps)
      results$future_max_temp[i] <- max(future_temps)
      
      # Calculate temperature increase
      results$temp_increase[i] <- results$future_mean_temp[i] - results$mean_temp[i]
      
      # Get thermal envelope from current climate
      max_threshold <- results$max_temp[i]
      min_threshold <- results$min_temp[i]
      
      # Calculate percentage exceeding thermal envelope
      above_max <- sum(future_temps > max_threshold)
      below_min <- sum(future_temps < min_threshold)
      total_points <- length(future_temps)
      
      results$exceed_max_pct[i] <- (above_max / total_points) * 100
      results$exceed_min_pct[i] <- (below_min / total_points) * 100
      results$exceed_envelope_pct[i] <- ((above_max + below_min) / total_points) * 100
    }
  }
  
  return(results)
}


# Set seed for reproducibility
set.seed(123)

# Run the simulation
today_temp <- temperature_gradient(max_temp=30, min_temp=0,latitude_range = c(0, 1000))
future_temp <- temperature_gradient(max_temp=33, min_temp=3,latitude_range = c(0, 1000))
segments <- generate_segments(1000,1,1000,latitude_range = c(0, 1000),distribution = "lognormal")
results <- analyze_temp_exceed(segments,today_temp,future_temp)



### Plot

today_temp$scenario <- "Today"
future_temp$scenario <- "Future"
combined_temp <- rbind(today_temp, future_temp)

combined_temp <- combined_temp %>%
  mutate(display_latitude = ((latitude - min(latitude)) / (max(latitude) - min(latitude))) * 100 - 50)

plot1 <- ggplot(combined_temp, aes(x = display_latitude, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Temperature (°C)",
       color = "Scenario")+
  scale_x_continuous(breaks = seq(-50, 50, by = 50)) 
  #theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

results$pct_domain <- results$size / 1000

# plot range size vs temperature change
plot2 <- ggplot(results, aes(x = pct_domain*100, y = temp_increase)) +
  geom_pointdensity(position = position_jitter(width = 0.01, height = 0.05),size=0.01)+
  #geom_hex()+
  scale_color_gradient(low='#e9f6fd',high="#2171b5") + 
  #geom_smooth(method = "loess", color = "black") +
  theme_classic() +
  labs(x = "Latitudinal range size (degrees)", 
       y = "Mean temperature \nincrease (°C)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20))+
  scale_y_continuous(limits=c(0,5))+
  theme(legend.position="none")

results <- results %>%
  arrange(pct_domain) %>%
  mutate(cum_species_pct = (row_number() / n()) * 100)

# plot range size vs habitat loss
plot3 <- ggplot(results, aes(x = pct_domain*100, y = exceed_max_pct)) +
  geom_pointdensity()+
  #geom_hex(bins=50)+
  scale_color_gradient(low='#e9f6fd',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
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
       y = "Habitat loss \n(% above today's max temp)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20))+
  theme(legend.position="none")

#plot range size distribution
plot4 <- ggplot(results,aes(x=pct_domain*100))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Latitudinal range size (degrees)", 
       y = "Frequency of species")+
  scale_x_continuous(limits = c(-5, 100),
                     breaks = seq(0, 100, 20))

plot_grid(plot1, plot2,plot3,plot4, labels = c("A", "B","C","D"), ncol = 4)










### climate displacement vs range size simulation

##functions

# Function to calculate latitudinal displacement of temperatures
calculate_latitudinal_displacement <- function(today_temp, future_temp) {
  # Sort both datasets by latitude
  today_temp <- today_temp[order(today_temp$latitude), ]
  future_temp <- future_temp[order(future_temp$latitude), ]
  
  # Create results dataframe based on today's points
  result <- data.frame(
    latitude = today_temp$latitude,
    today_temp = today_temp$temperature,
    future_temp = future_temp$temperature,
    future_lat = NA_real_,
    future_temp_analog = NA_real_,
    displacement = NA_real_
  )
  
  # For each today's temperature point
  for (i in 1:nrow(result)) {
    current_lat <- result$latitude[i]
    current_temp <- result$today_temp[i]
    
    # Skip if temperature is NA
    if (is.na(current_temp)) {
      next
    }
    
    # Find all points in future data that have approximately this temperature
    # Use a small tolerance to find close matches
    tolerance <- 0.1  # Adjust based on your temperature precision
    matches <- which(abs(future_temp$temperature - current_temp) < tolerance)
    
    # If no matches within tolerance, check if it's within the range
    if (length(matches) == 0) {
      min_future_temp <- min(future_temp$temperature, na.rm = TRUE)
      max_future_temp <- max(future_temp$temperature, na.rm = TRUE)
      
      # Skip if outside temperature range
      if (current_temp < min_future_temp || current_temp > max_future_temp) {
        next
      }
      
      # Find closest temperature points
      temp_diffs <- abs(future_temp$temperature - current_temp)
      sorted_indices <- order(temp_diffs)
      
      # Get the two closest points
      idx1 <- sorted_indices[1]
      idx2 <- sorted_indices[2]
      
      # Points should be on opposite sides of the target temperature
      if (sign(future_temp$temperature[idx1] - current_temp) == 
          sign(future_temp$temperature[idx2] - current_temp)) {
        # Find another point on the opposite side
        for (j in 3:length(sorted_indices)) {
          if (sign(future_temp$temperature[sorted_indices[j]] - current_temp) != 
              sign(future_temp$temperature[idx1] - current_temp)) {
            idx2 <- sorted_indices[j]
            break
          }
        }
      }
      
      # Interpolate between these points
      x1 <- future_temp$temperature[idx1]
      x2 <- future_temp$temperature[idx2]
      y1 <- future_temp$latitude[idx1]
      y2 <- future_temp$latitude[idx2]
      
      # Linear interpolation
      if (abs(x2 - x1) < 1e-6) {
        result$future_lat[i] <- (y1 + y2) / 2
      } else {
        result$future_lat[i] <- y1 + (y2 - y1) * (current_temp - x1) / (x2 - x1)
      }
      result$future_temp[i] <- current_temp
    } else {
      # If we have exact matches, find the one closest in latitude
      lat_diffs <- abs(future_temp$latitude[matches] - current_lat)
      closest_match <- matches[which.min(lat_diffs)]
      
      result$future_lat[i] <- future_temp$latitude[closest_match]
      result$future_temp_analog[i] <- future_temp$temperature[closest_match]
    }
    
    # Calculate displacement
    result$displacement[i] <- result$future_lat[i] - result$latitude[i]
  }
  
  # Create display-ready values
  min_lat <- min(c(today_temp$latitude, future_temp$latitude), na.rm = TRUE)
  max_lat <- max(c(today_temp$latitude, future_temp$latitude), na.rm = TRUE)
  lat_range <- max_lat - min_lat
  
  # Scale to -50 to 50 for display (changed from -90 to 90)
  result$display_latitude <- -50 + 100 * (result$latitude - min_lat) / lat_range
  
  # Scale displacement proportionally
  scale_factor <- 100 / lat_range  # Changed from 180 to 100
  result$display_displacement <- abs(result$displacement * scale_factor)
  
  return(result)
}

# Create a function to identify critical range size at each latitude
find_critical_range_size <- function(displacement_data, segments_data) {
  # Clone segments data
  segments_display <- segments_data
  
  # Calculate latitude range from displacement data
  lat_min <- min(displacement_data$latitude, na.rm=TRUE)
  lat_max <- max(displacement_data$latitude, na.rm=TRUE)
  lat_range <- lat_max - lat_min
  
  # Keep original segment sizes
  segments_display$size <- segments_display$end - segments_display$start
  
  # Calculate center point for each segment
  segments_display$center_lat <- (segments_display$start + segments_display$end) / 2
  
  # Calculate 1/4 and 3/4 points
  segments_display$quarter_lat <- (segments_display$start + segments_display$center_lat) / 2
  segments_display$three_quarter_lat <- (segments_display$center_lat + segments_display$end) / 2
  
  # Initialize results columns
  segments_display$required_displacement <- numeric(nrow(segments_display))
  segments_display$min_required_displacement <- numeric(nrow(segments_display))
  segments_display$no_analog <- logical(nrow(segments_display))
  
  for (i in 1:nrow(segments_display)) {
    # Get all five points for this segment
    points <- c(
      segments_display$start[i],
      segments_display$quarter_lat[i],
      segments_display$center_lat[i],
      segments_display$three_quarter_lat[i],
      segments_display$end[i]
    )
    
    # Store displacements for each point
    displacements <- numeric(length(points))
    has_analog <- logical(length(points))
    
    # Find displacement for each point
    for (j in 1:length(points)) {
      point <- points[j]
      closest_idx <- which.min(abs(displacement_data$latitude - point))
      
      # Check if displacement is NA
      if (is.na(displacement_data$displacement[closest_idx])) {
        displacements[j] <- NA
        has_analog[j] <- FALSE
      } else {
        displacements[j] <- abs(displacement_data$displacement[closest_idx])
        has_analog[j] <- TRUE
      }
    }
    
    # Store the center point displacement for reference (original behavior)
    segments_display$required_displacement[i] <- displacements[3]  # Center point
    
    # Find the minimum displacement among points that have analogs
    valid_displacements <- displacements[has_analog]
    if (length(valid_displacements) > 0) {
      segments_display$min_required_displacement[i] <- min(valid_displacements, na.rm=TRUE)
      segments_display$no_analog[i] <- FALSE
    } else {
      segments_display$min_required_displacement[i] <- NA
      segments_display$no_analog[i] <- TRUE
    }
  }
  
  # Determine if the segment will completely lose habitat using the minimum displacement
  # If min displacement is less than segment size, at least some portion can adapt
  segments_display$complete_habitat_loss <- segments_display$size < segments_display$min_required_displacement
  
  # Areas with no analog should be marked as complete habitat loss
  segments_display$complete_habitat_loss[segments_display$no_analog] <- TRUE
  
  # Convert to display units at the end - changed from -90/90 to -50/50
  segments_display$display_start <- ((segments_display$start - 0) / lat_range) * 100 - 50
  segments_display$display_end <- ((segments_display$end - 0) / lat_range) * 100 - 50
  segments_display$display_size <- segments_display$display_end - segments_display$display_start
  segments_display$display_center_lat <- ((segments_display$center_lat - 0) / lat_range) * 100 - 50
  
  # Scale displacement to display units - changed scale factor from 180 to 100
  scale_factor <- 100 / lat_range
  segments_display$display_displacement <- segments_display$required_displacement * scale_factor
  segments_display$display_min_displacement <- segments_display$min_required_displacement * scale_factor
  
  return(segments_display)
}


## run

# Calculate latitudinal displacement
displacement <- calculate_latitudinal_displacement(today_temp, future_temp)

# Apply this to your segments
analyzed_segments <- find_critical_range_size(displacement, segments)

# Calculate the percentage of segments that will face complete habitat loss
percent_loss <- sum(analyzed_segments$complete_habitat_loss) / nrow(analyzed_segments) * 100
analyzed_segments$pct_domain <- analyzed_segments$size / 1000

combined_segment_results <- merge(results, analyzed_segments, by = "segment_id", all.x = TRUE)

segments_lost <- subset(analyzed_segments, complete_habitat_loss == TRUE)

# Calculate statistics about these segments
critical_stats <- data.frame(
  mean_size = mean(segments_lost$size, na.rm = TRUE),
  median_size = median(segments_lost$size, na.rm = TRUE),
  max_size = max(segments_lost$size, na.rm = TRUE)
)

# Use the maximum size as the critical threshold
# This represents the largest range that still experiences complete habitat loss
critical_range_size <- median(segments_lost$size, na.rm = TRUE)
critical_range_pct <- (critical_range_size / max(analyzed_segments$size, na.rm = TRUE)) * 100

## plot

# Plot displacement vs. latitude
plot_displacement <- ggplot(displacement, aes(x = display_latitude, y = display_displacement)) +
  geom_line(size = 1) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Displacement (degrees latitude)",
       title = "Required latitudinal shift to maintain current temperature")+
  scale_x_continuous(limits = c(-50, 50))+
  scale_y_continuous(limits=c(0,30))

# Plot showing critical threshold
plot_critical <- ggplot() +
  # Standard points with displacement values
  geom_point(data = subset(analyzed_segments, !no_analog), 
             aes(x = display_size, y = display_displacement, 
                 color = complete_habitat_loss), 
             alpha = 0.5) +
  # No analog points jittered near x-axis
  geom_point(data = subset(analyzed_segments, no_analog),
             aes(x = display_size, y = 0), 
             color = "purple", 
             position = position_jitter(height = 0.5), 
             alpha = 0.5) +
  # Critical threshold line
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  # Colors
  scale_color_manual(values = c("blue", "red"), name = "Complete habitat loss") +
  theme_classic() +
  labs(x = "Range size (degrees latitude)", 
       y = "Required displacement (degrees latitude)") 

#plot3 with critical range size shading
min_critical_range <- min(displacement$display_displacement, na.rm = TRUE)
max_critical_range <- max(displacement$display_displacement, na.rm = TRUE)

plot3_crit <- plot3 +  annotate("rect", 
                                xmin = min_critical_range, 
                                xmax = max_critical_range, 
                                ymin = -Inf, 
                                ymax = Inf, 
                                fill = "red", 
                                alpha = 0.2) 


# Plot showing habitat loss by latitude
# Calculate percentage of segments at risk
future_temp6 <- temperature_gradient(max_temp=36, min_temp=6,latitude_range = c(0, 1000))
displacement6 <- calculate_latitudinal_displacement(today_temp, future_temp6)
future_temp9 <- temperature_gradient(max_temp=39, min_temp=9,latitude_range = c(0, 1000))
displacement9 <- calculate_latitudinal_displacement(today_temp, future_temp9)

plot_by_latitude <- ggplot(analyzed_segments, aes(x = display_center_lat, y = display_size)) +
  # Regular points
  geom_point(aes(color = complete_habitat_loss), size = 0.5,alpha = 0.5) +
  # No analog points 
  geom_point(data = subset(analyzed_segments, no_analog), 
             aes(x = display_center_lat, y = display_size),
             color = "red", size = 0.5, alpha = 0.5) +
  # required latitudinal shift
  geom_line(data = displacement, aes(x = display_latitude, y = abs(display_displacement)), 
            color = "black", size = 1.5) +
  geom_line(data = displacement6, aes(x = display_latitude, y = abs(display_displacement)), 
            color = "darkgrey", size = 1.5) +
  geom_line(data = displacement9, aes(x = display_latitude, y = abs(display_displacement)), 
            color = "lightgrey", size = 1.5) +
  scale_color_manual(values = c("blue", "red"), name = "Complete \nhabitat loss") +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") 


# Create habitat vulnerability plot with threshold line
habitat_vulnerability_plot <- ggplot(combined_segment_results, aes(x = pct_domain.x * 100, y = exceed_envelope_pct)) +
  # Regular points colored by habitat loss (blue = FALSE, red = TRUE)
  geom_point(aes(color = complete_habitat_loss)) +
  
  # Add vertical line at critical threshold
  geom_vline(xintercept = critical_range_pct, 
             color = "red", linewidth = 1, linetype = "dashed") +
  
  # Highlight no analog points with purple triangles
  geom_point(data = subset(combined_segment_results, no_analog == TRUE), 
             aes(x = pct_domain.x * 100, y = exceed_envelope_pct),
             color = "purple", shape = 17, size = 3,alpha=0.1) +
  
  # Set colors - blue for false, red for true
  scale_color_manual(values = c("blue", "red"), 
                     name = "Range too small\nfor adaptation") +
  
  # Labels
  labs(x = "Latitudinal range (% of domain)",
       y = "Habitat loss (% above today's max temp)") +
  theme_classic() 



# bar chart: proportion of species in critical range size
percent_loss <- sum(analyzed_segments$complete_habitat_loss, na.rm = TRUE) / nrow(analyzed_segments) * 100

points_under_line <- function(data, displacement_data) {
  count <- 0
  # Only count rows that have complete data (no NAs)
  complete_data <- data[complete.cases(data$display_center_lat, data$display_size), ]
  total <- nrow(complete_data)
  
  for(i in 1:total) {
    # Find the closest latitude in displacement data
    lat <- complete_data$display_center_lat[i]
    closest_idx <- which.min(abs(displacement_data$display_latitude - lat))
    
    # Make sure we have valid threshold value
    if(length(closest_idx) > 0 && !is.na(displacement_data$display_displacement[closest_idx])) {
      threshold <- displacement_data$display_displacement[closest_idx]
      
      # Check if point is under the line
      if(!is.na(complete_data$display_size[i]) && !is.na(threshold) && 
         complete_data$display_size[i] < abs(threshold)) {
        count <- count + 1
      }
    }
  }
  
  return(count/total * 100)
}
pct_under_line1 <- points_under_line(analyzed_segments, displacement)
pct_under_line2 <- points_under_line(analyzed_segments, displacement6)
pct_under_line3 <- points_under_line(analyzed_segments, displacement9)

# Create a data frame for plotting
bar_data <- data.frame(
  Line = c("3C", "6C", "9C"),
  Percentage = c(pct_under_line1, pct_under_line2, pct_under_line3)
)

# Create bar chart
barE <- ggplot(bar_data, aes(x = Line, y = Percentage)) +
  geom_bar(stat = "identity", fill = c("black", "darkgrey", "lightgrey")) +
  theme_classic() +
  labs(x = "Displacement line", 
       y = "Critical range sizes (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(limits=c(0,100))

plot_grid(plot1, plot2,plot3_crit,plot4,plot_by_latitude,barE,
          labels = c("A", "B","C","D","E","F"), ncol = 2,nrow=3)






















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
generate_polygons <- function(num_polygons, grid_size, min_size, max_size = grid_size * 0.9) {
  # Create lists to store polygon data
  all_polygons <- list()
  ids <- numeric(num_polygons)
  
  # Center longitude for all polygons
  center_x <- grid_size / 2
  
  # Distribute polygons across size categories
  num_very_small <- floor(num_polygons * 0.3)  # 30% are very small (<5%)
  num_small <- floor(num_polygons * 0.25)      # 25% are small (5-20%)
  num_medium <- floor(num_polygons * 0.2)     # 15% are medium (20-50%)
  num_large <- num_polygons - num_very_small - num_small - num_medium  # 10% are large (50-90%)
  
  polygon_count <- 0
  
  # 1. Generate very small polygons (0-5% of domain)
  for (i in 1:num_very_small) {
    polygon_count <- polygon_count + 1
    
    # Size between min_size and 5% of domain
    max_radius <- grid_size * 0.025  # 5% diameter = 2.5% radius
    radius <- min_size + runif(1) * (max_radius - min_size)
    
    # Random latitude within safe bounds
    min_y <- radius * 1.2
    max_y <- grid_size - (radius * 1.2)
    center_y <- runif(1, min_y, max_y)
    
    # Generate polygon
    num_vertices <- sample(5:8, 1)
    angles <- sort(runif(num_vertices, 0, 2*pi))
    radii <- radius * runif(num_vertices, 0.7, 1.0)
    
    x_coords <- center_x + radii * cos(angles)
    y_coords <- center_y + radii * sin(angles)
    
    # Close polygon
    x_coords <- c(x_coords, x_coords[1])
    y_coords <- c(y_coords, y_coords[1])
    
    coords <- cbind(x_coords, y_coords)
    all_polygons[[polygon_count]] <- list(coords)
    ids[polygon_count] <- polygon_count
  }
  
  # 2. Generate small polygons (5-20% of domain)
  for (i in 1:num_small) {
    polygon_count <- polygon_count + 1
    
    # Size between 5% and 20% of domain
    min_radius <- grid_size * 0.025  # 5% diameter = 2.5% radius
    max_radius <- grid_size * 0.1    # 20% diameter = 10% radius
    radius <- min_radius + runif(1) * (max_radius - min_radius)
    
    # Random latitude within safe bounds
    min_y <- radius * 1.2
    max_y <- grid_size - (radius * 1.2)
    
    # Skip if invalid bounds
    if (min_y >= max_y) {
      radius <- (max_y - 10) / 1.2
    }
    
    center_y <- runif(1, min_y, max_y)
    
    # Generate polygon
    num_vertices <- sample(6:10, 1)
    angles <- sort(runif(num_vertices, 0, 2*pi))
    radii <- radius * runif(num_vertices, 0.7, 1.0)
    
    x_coords <- center_x + radii * cos(angles)
    y_coords <- center_y + radii * sin(angles)
    
    # Close polygon
    x_coords <- c(x_coords, x_coords[1])
    y_coords <- c(y_coords, y_coords[1])
    
    coords <- cbind(x_coords, y_coords)
    all_polygons[[polygon_count]] <- list(coords)
    ids[polygon_count] <- polygon_count
  }
  
  # 3. Generate medium polygons (20-50% of domain)
  for (i in 1:num_medium) {
    polygon_count <- polygon_count + 1
    
    # Size between 20% and 50% of domain
    min_radius <- grid_size * 0.1    # 20% diameter = 10% radius
    max_radius <- grid_size * 0.25   # 50% diameter = 25% radius
    radius <- min_radius + runif(1) * (max_radius - min_radius)
    
    # More restricted placement for medium polygons
    max_offset <- (grid_size / 2) - radius
    
    # If radius too large, adjust
    if (max_offset < 0) {
      radius <- (grid_size / 2) * 0.9
      max_offset <- (grid_size / 2) - radius
    }
    
    # Center near middle of grid
    center_y <- (grid_size / 2) + runif(1, -max_offset, max_offset)
    
    # Generate polygon
    num_vertices <- sample(8:12, 1)
    angles <- sort(runif(num_vertices, 0, 2*pi))
    radii <- radius * runif(num_vertices, 0.8, 1.0)
    
    x_coords <- center_x + radii * cos(angles)
    y_coords <- center_y + radii * sin(angles)
    
    # Close polygon
    x_coords <- c(x_coords, x_coords[1])
    y_coords <- c(y_coords, y_coords[1])
    
    # Clip to grid if needed
    x_coords <- pmin(pmax(x_coords, 0), grid_size)
    y_coords <- pmin(pmax(y_coords, 0), grid_size)
    
    coords <- cbind(x_coords, y_coords)
    all_polygons[[polygon_count]] <- list(coords)
    ids[polygon_count] <- polygon_count
  }
  
  # 4. Generate large polygons (50-90% of domain)
  for (i in 1:num_large) {
    polygon_count <- polygon_count + 1
    
    # Size between 50% and 90% of domain
    # Spread them out evenly across this range
    size_pct <- 0.5 + (0.4 * i / num_large)
    radius <- (grid_size / 2) * size_pct
    
    # Large polygons must be very close to center
    max_offset <- max(0, (grid_size / 2) - radius)
    center_y <- (grid_size / 2) + runif(1, -max_offset, max_offset)
    
    # More vertices for smoother large polygons
    num_vertices <- 16
    angles <- sort(runif(num_vertices, 0, 2*pi))
    
    # Less variation for large polygons
    variation <- 0.95
    radii <- radius * runif(num_vertices, variation, 1.0)
    
    x_coords <- center_x + radii * cos(angles)
    y_coords <- center_y + radii * sin(angles)
    
    # Close polygon
    x_coords <- c(x_coords, x_coords[1])
    y_coords <- c(y_coords, y_coords[1])
    
    # Clip to grid boundaries
    x_coords <- pmin(pmax(x_coords, 0), grid_size)
    y_coords <- pmin(pmax(y_coords, 0), grid_size)
    
    coords <- cbind(x_coords, y_coords)
    all_polygons[[polygon_count]] <- list(coords)
    ids[polygon_count] <- polygon_count
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
    segment_spacing = 0.1  # Add spacing between segments (as a fraction of segment width)
) {
  # Get latitude range from data
  latitude_min <- min(latitude_sample$y)
  latitude_max <- max(latitude_sample$y)
  
  # Split the data by scenario
  today_data <- latitude_sample %>% filter(scenario == "Today")
  future_data <- latitude_sample %>% filter(scenario == "Future")
  
  # Initialize output dataframe for segments
  all_segments <- data.frame()
  
  # Generate segments for each row
  for (row in 1:length(segments_per_row)) {
    num_segments <- segments_per_row[row]
    
    # Position for this row
    y_position <- bottom_margin - ((row - 1) * y_spacing)
    
    if (num_segments == 1) {
      # For the single segment row, don't apply spacing - use full width
      adj_range_min <- latitude_min
      adj_range_max <- latitude_max
    } else {
      # For other rows, calculate segment width with spacing as before
      segment_width <- (latitude_max - latitude_min) / num_segments
      effective_width <- segment_width * (1 - segment_spacing)
      
      # Create breakpoints dividing the latitude range into equal segments
      breakpoints <- seq(latitude_min, latitude_max, length.out = num_segments + 1)
    }
    
    # Create segments
    if (num_segments == 1) {
      # For the single segment row, process the whole range at once
      range_min <- latitude_min
      range_max <- latitude_max
      
      # For envelope determination, use dense sampling across the range
      sample_points <- 100
      latitude_points <- seq(range_min, range_max, length.out = sample_points)
      
      # Get temperatures at these points
      today_temps <- approx(today_data$y, today_data$temperature, latitude_points)$y
      future_temps <- approx(future_data$y, future_data$temperature, latitude_points)$y
      
      # Find min/max of today's temperatures for this range
      min_envelope <- min(today_temps, na.rm = TRUE)
      max_envelope <- max(today_temps, na.rm = TRUE)
      
      # Check if future temperatures are within envelope
      within_envelope <- future_temps >= min_envelope & future_temps <= max_envelope
      
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
      
      # For a single segment, use the exact latitudes without adjusting
      vis_transitions <- transition_points
      
      # Create segments
      for (j in 1:(length(vis_transitions)-1)) {
        # Create segment from current transition to next
        all_segments <- rbind(all_segments, data.frame(
          row_id = row,
          segment_id = paste(row, 1, j, sep = "-"),
          y_start = vis_transitions[j],
          y_end = vis_transitions[j+1],
          y_position = y_position,
          within_envelope = if(j <= length(rle_result$values)) rle_result$values[j] else !rle_result$values[length(rle_result$values)]
        ))
      }
    } else {
      # For multi-segment rows, process each segment separately
      for (i in 1:num_segments) {
        range_min <- breakpoints[i]
        range_max <- breakpoints[i + 1]
        
        # Apply spacing by shrinking segment width
        range_center <- (range_min + range_max) / 2
        adj_range_min <- range_center - (effective_width / 2)
        adj_range_max <- range_center + (effective_width / 2)
        
        # For envelope determination, use dense sampling for accurate checking
        sample_points <- 100
        latitude_points <- seq(range_min, range_max, length.out = sample_points)
        
        # Get temperatures at these points
        today_temps <- approx(today_data$y, today_data$temperature, latitude_points)$y
        future_temps <- approx(future_data$y, future_data$temperature, latitude_points)$y
        
        # Find min/max of today's temperatures for this range
        min_envelope <- min(today_temps, na.rm = TRUE)
        max_envelope <- max(today_temps, na.rm = TRUE)
        
        # Check if future temperatures are within envelope
        within_envelope <- future_temps >= min_envelope & future_temps <= max_envelope
        
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
        
        # Map to adjusted (visualization) coordinates
        vis_transitions <- adj_range_min + (adj_range_max - adj_range_min) * 
          (transition_points - range_min) / (range_max - range_min)
        
        # Create segments
        for (j in 1:(length(vis_transitions)-1)) {
          # Create segment from current transition to next
          all_segments <- rbind(all_segments, data.frame(
            row_id = row,
            segment_id = paste(row, i, j, sep = "-"),
            y_start = vis_transitions[j],
            y_end = vis_transitions[j+1],
            y_position = y_position,
            within_envelope = if(j <= length(rle_result$values)) rle_result$values[j] else !rle_result$values[length(rle_result$values)]
          ))
        }
      }
    }
  }
  
  return(all_segments)
}


create_latitude_segments <- function(
    latitude_sample,  # The latitude_sample data frame from your existing code
    segments_per_row = c(12, 6, 3, 1),  # Number of segments for each row
    y_spacing = 5,    # Vertical spacing between rows
    bottom_margin = -10, # Position of the bottom-most row
    segment_spacing = 0.1  # Add spacing between segments (as a fraction of segment width)
) {
  # Get latitude range from data
  latitude_min <- min(latitude_sample$y)
  latitude_max <- max(latitude_sample$y)
  
  # Split the data by scenario
  today_data <- latitude_sample %>% filter(scenario == "Today")
  future_data <- latitude_sample %>% filter(scenario == "Future")
  
  # Initialize output dataframe for segments
  all_segments <- data.frame()
  
  # Initialize dataframe for row statistics
  row_stats <- data.frame(
    row_id = integer(),
    avg_temp_increase = numeric(),
    pct_exceeding_envelope = numeric()
  )
  
  # Generate segments for each row
  for (row in 1:length(segments_per_row)) {
    num_segments <- segments_per_row[row]
    
    # Position for this row
    y_position <- bottom_margin - ((row - 1) * y_spacing)
    
    # Track statistics for this row
    row_segments <- data.frame()
    total_width <- 0
    weighted_temp_increase <- 0
    total_exceeding_width <- 0
    
    # MODIFIED: Special handling for the last row (1 segment)
    if (num_segments == 1) {
      # For the single segment row, don't apply spacing - use full width
      adj_range_min <- latitude_min
      adj_range_max <- latitude_max
    } else {
      # For other rows, calculate segment width with spacing as before
      segment_width <- (latitude_max - latitude_min) / num_segments
      effective_width <- segment_width * (1 - segment_spacing)
      
      # Create breakpoints dividing the latitude range into equal segments
      breakpoints <- seq(latitude_min, latitude_max, length.out = num_segments + 1)
    }
    
    # Create segments
    if (num_segments == 1) {
      # For the single segment row, process the whole range at once
      range_min <- latitude_min
      range_max <- latitude_max
      
      # For envelope determination, use dense sampling across the range
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
      
      # For a single segment, use the exact latitudes without adjusting
      vis_transitions <- transition_points
      
      # Create segments
      for (j in 1:(length(vis_transitions)-1)) {
        segment_width <- vis_transitions[j+1] - vis_transitions[j]
        is_within <- if(j <= length(rle_result$values)) rle_result$values[j] else !rle_result$values[length(rle_result$values)]
        
        # Create segment from current transition to next
        segment_data <- data.frame(
          row_id = row,
          segment_id = paste(row, 1, j, sep = "-"),
          y_start = vis_transitions[j],
          y_end = vis_transitions[j+1],
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
        if (!is_within) {
          total_exceeding_width <- total_exceeding_width + segment_width
        }
      }
    } else {
      # For multi-segment rows, process each segment separately
      for (i in 1:num_segments) {
        range_min <- breakpoints[i]
        range_max <- breakpoints[i + 1]
        
        # Apply spacing by shrinking segment width
        range_center <- (range_min + range_max) / 2
        adj_range_min <- range_center - (effective_width / 2)
        adj_range_max <- range_center + (effective_width / 2)
        
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
        
        # Map to adjusted (visualization) coordinates
        vis_transitions <- adj_range_min + (adj_range_max - adj_range_min) * 
          (transition_points - range_min) / (range_max - range_min)
        
        # Create segments
        for (j in 1:(length(vis_transitions)-1)) {
          segment_width <- vis_transitions[j+1] - vis_transitions[j]
          is_within <- if(j <= length(rle_result$values)) rle_result$values[j] else !rle_result$values[length(rle_result$values)]
          
          # Create segment from current transition to next
          segment_data <- data.frame(
            row_id = row,
            segment_id = paste(row, i, j, sep = "-"),
            y_start = vis_transitions[j],
            y_end = vis_transitions[j+1],
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
          if (!is_within) {
            total_exceeding_width <- total_exceeding_width + segment_width
          }
        }
      }
    }
    
    # Calculate row statistics
    avg_temp_increase <- weighted_temp_increase / total_width
    pct_exceeding_envelope <- (total_exceeding_width / total_width) * 100
    
    # Add row statistics
    row_stats <- rbind(row_stats, data.frame(
      row_id = row,
      avg_temp_increase = avg_temp_increase,
      pct_exceeding_envelope = pct_exceeding_envelope
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

# Create future temperature gradient
future_temp_raster <- create_temperature_gradient(grid_size, max_temp=35, min_temp=5)

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
plot1 <- ggplot(latitude_sample, aes(x = y, y = temperature, color = scenario)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Today" = "blue", "Future" = "red")) +
  theme_classic() +
  labs(x = "Latitude", 
       y = "Temperature (°C)",
       color = "Scenario")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


segment_result <- create_latitude_segments(latitude_sample)
segments_plot <- segment_result$segments
row_stats <- segment_result$stats

y_position_map <- c(
  "-10" = -4,  # Current -10 -> new -2.5
  "-15" = -5,    # Current -15 -> new -5
  "-20" = -6,  # Current -20 -> new -7.5
  "-25" = -7    # Current -25 -> new -10
)

# Apply the mapping to create a new y_position column
segments_plot$y_position_new <- y_position_map[as.character(segments_plot$y_position)]

# Format the statistics for display
stats_text <- row_stats %>%
  mutate(
    label = sprintf("%.1f°C, %.1f%%",
                    avg_temp_increase, 
                    pct_exceeding_envelope))


plot1 + 
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






plot2_2d <- ggplot(results_2d, aes(x = area/ (grid_size * grid_size), y = temp_increase)) +
  #geom_point(alpha = 0.6, color = "darkred") +
  geom_pointdensity()+
  #geom_hex()+
  scale_color_gradient(low='#eff3ff',high="#2171b5") + 
  #geom_smooth(method = "loess", color = "black") +
  theme_classic() +
  labs(x = "% of domain", 
       y = "Mean temperature increase (°C)")+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100),
                     breaks = seq(0, 1, 0.2))+
  theme(legend.position="none")


# plot range size vs habitat loss
plot3_2d <- ggplot(results_2d, aes(x = area / (grid_size * grid_size), y = exceed_max_pct)) +
  #geom_pointdensity()+
  geom_hex()+
  scale_fill_gradient(low='#eff3ff',high="#2171b5") + #scale_color for point, scale_fill for hex
  #scale_fill_viridis_c()+
  stat_summary_bin(
    fun = "mean", 
    geom = "line",
    color = "black",
    linewidth = 1,
    bins = 10  
  ) +
  #geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  theme_classic() +
  labs(x = "Latitudinal range", 
       y = "Habitat loss \n(% above today's max temp)")+
  scale_x_continuous(limits = c(0, 1), 
                     labels = scales::percent_format(scale = 100))+
  theme(legend.position="none")

#plot range size distribution
plot4_2d <- ggplot(results_2d,aes(x=area / (grid_size * grid_size)))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Latitudinal range", 
       y = "Count")+
  scale_x_continuous(limits = c(-0.05, 1), 
                     labels = scales::percent_format(scale = 100))

plot_grid(plot1, plot2_2d,plot3_2d,plot4_2d, labels = c("A", "B","C","D"), ncol = 4)


plot_grid(plot1, plot2,plot3, labels = c("A", "B","C"), ncol = 3)


