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
  results_sim <- data.frame(
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
    results_sim$mean_temp[i] <- mean(segment_temps)
    results_sim$min_temp[i] <- min(segment_temps)
    results_sim$max_temp[i] <- max(segment_temps)
  }
  
  # If future gradient is provided, calculate additional metrics
  if (!is.null(future_gradient)) {
    # Add future temperature columns
    results_sim$future_mean_temp <- numeric(nrow(segments))
    results_sim$future_min_temp <- numeric(nrow(segments))
    results_sim$future_max_temp <- numeric(nrow(segments))
    results_sim$temp_increase <- numeric(nrow(segments))
    results_sim$exceed_max_pct <- numeric(nrow(segments))
    results_sim$exceed_min_pct <- numeric(nrow(segments))
    results_sim$exceed_envelope_pct <- numeric(nrow(segments))
    
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
      results_sim$future_mean_temp[i] <- mean(future_temps)
      results_sim$future_min_temp[i] <- min(future_temps)
      results_sim$future_max_temp[i] <- max(future_temps)
      
      # Calculate temperature increase
      results_sim$temp_increase[i] <- results_sim$future_mean_temp[i] - results_sim$mean_temp[i]
      
      # Get thermal envelope from current climate
      max_threshold <- results_sim$max_temp[i]
      min_threshold <- results_sim$min_temp[i]
      
      # Calculate percentage exceeding thermal envelope
      above_max <- sum(future_temps > max_threshold)
      below_min <- sum(future_temps < min_threshold)
      total_points <- length(future_temps)
      
      results_sim$exceed_max_pct[i] <- (above_max / total_points) * 100
      results_sim$exceed_min_pct[i] <- (below_min / total_points) * 100
      results_sim$exceed_envelope_pct[i] <- ((above_max + below_min) / total_points) * 100
    }
  }
  
  return(results_sim)
}


# Set seed for reproducibility
set.seed(123)

# Run the simulation
today_temp <- temperature_gradient(max_temp=30, min_temp=0,latitude_range = c(0, 1000))
future_temp <- temperature_gradient(max_temp=33, min_temp=3,latitude_range = c(0, 1000))
segments <- generate_segments(1000,1,1000,latitude_range = c(0, 1000),distribution = "lognormal")
results_sim <- analyze_temp_exceed(segments,today_temp,future_temp)



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
  scale_x_continuous(breaks = seq(-50, 50, by = 50)) +
  theme(legend.position = c(0.5, 0.1),  # x,y coordinates (0.5 is center, 0.1 is near bottom)
        legend.justification = c(0.5, 0), # anchor point for the legend
        legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5),
        legend.margin = margin(6, 6, 6, 6))
  #theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

results_sim$pct_domain <- results_sim$size / 1000

# plot range size vs temperature change
plot2 <- ggplot(results_sim, aes(x = pct_domain*100, y = temp_increase)) +
  geom_pointdensity(position = position_jitter(width = 0.01, height = 0.05),size=0.01)+
  scale_color_gradient(low='#e9f6fd',high="#2171b5") + 
  theme_classic() +
  labs(x = "Range size as\nproportion of domain", 
       y = "Temperature \nexposure (°C)")+
  scale_x_continuous(limits = c(0, 100),
                     expand=c(0,0),
                     breaks = seq(0, 100, 20))+
  scale_y_continuous(limits=c(0,5),
                     expand=c(0,0),)+
  theme(legend.position="none")

results_sim <- results_sim %>%
  arrange(pct_domain) %>%
  mutate(cum_species_pct = (row_number() / n()) * 100)



#plot range size distribution
plot4 <- ggplot(results_sim,aes(x=pct_domain*100))+
  geom_histogram(fill = "lightgray", color = "black")+
  theme_classic()+
  labs(x = "Range size as\nproportion of domain", 
       y = "Frequency \nof species")+
  scale_x_continuous(limits = c(-5, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  scale_y_continuous(limits=c(0,300),
                     expand=c(0,0))











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
  
  # Keep original calculations for backward compatibility
  segments_display$center_lat <- (segments_display$start + segments_display$end) / 2
  segments_display$quarter_lat <- (segments_display$start + segments_display$center_lat) / 2
  segments_display$three_quarter_lat <- (segments_display$center_lat + segments_display$end) / 2
  
  # Initialize results columns
  segments_display$required_displacement <- numeric(nrow(segments_display))
  segments_display$min_required_displacement <- numeric(nrow(segments_display))
  segments_display$no_analog <- logical(nrow(segments_display))
  segments_display$points_sampled <- numeric(nrow(segments_display))
  segments_display$points_with_analog_anywhere <- numeric(nrow(segments_display))
  segments_display$points_with_analog_in_range <- numeric(nrow(segments_display))
  segments_display$no_analog_within_range_pct <- numeric(nrow(segments_display))
  
  for (i in 1:nrow(segments_display)) {
    # Find all data points within this segment instead of just 5 fixed points
    segment_indices <- which(displacement_data$latitude >= segments_display$start[i] & 
                               displacement_data$latitude <= segments_display$end[i])
    
    # If no points found in this range, fall back to the original 5-point method
    if (length(segment_indices) == 0) {
      points <- c(
        segments_display$start[i],
        segments_display$quarter_lat[i],
        segments_display$center_lat[i],
        segments_display$three_quarter_lat[i],
        segments_display$end[i]
      )
    } else {
      # Use all data points within the segment
      points <- displacement_data$latitude[segment_indices]
    }
    
    # Record how many points we're sampling
    segments_display$points_sampled[i] <- length(points)
    
    # Store displacements for each point
    displacements <- numeric(length(points))
    has_analog <- logical(length(points))
    within_range <- logical(length(points))
    
    # Find displacement for each point
    for (j in 1:length(points)) {
      point <- points[j]
      
      # Find the two closest points in displacement_data for interpolation
      dists <- abs(displacement_data$latitude - point)
      closest_idx <- which.min(dists)
      
      # If we're exactly at a data point, just use that value
      if (dists[closest_idx] < 1e-10) {
        if (is.na(displacement_data$displacement[closest_idx])) {
          displacements[j] <- NA
          has_analog[j] <- FALSE
          within_range[j] <- FALSE
        } else {
          displacements[j] <- displacement_data$displacement[closest_idx]
          has_analog[j] <- TRUE
          
          # Check if analog is within range
          future_lat <- point + displacements[j]
          within_range[j] <- (future_lat >= segments_display$start[i] && 
                                future_lat <= segments_display$end[i])
        }
      } else {
        # Get the two points surrounding our target for interpolation
        left_idx <- max(which(displacement_data$latitude <= point), na.rm = TRUE)
        right_idx <- min(which(displacement_data$latitude >= point), na.rm = TRUE)
        
        # Handle edge cases
        if (length(left_idx) == 0) left_idx <- right_idx
        if (length(right_idx) == 0) right_idx <- left_idx
        
        # Check for NA values
        if (is.na(displacement_data$displacement[left_idx]) || 
            is.na(displacement_data$displacement[right_idx])) {
          displacements[j] <- NA
          has_analog[j] <- FALSE
          within_range[j] <- FALSE
        } else {
          # Interpolate
          left_lat <- displacement_data$latitude[left_idx]
          right_lat <- displacement_data$latitude[right_idx]
          left_disp <- displacement_data$displacement[left_idx]
          right_disp <- displacement_data$displacement[right_idx]
          
          # If same point, just use that value
          if (left_idx == right_idx) {
            displacements[j] <- left_disp
          } else {
            # Linear interpolation
            weight <- (point - left_lat) / (right_lat - left_lat)
            displacements[j] <- left_disp + weight * (right_disp - left_disp)
          }
          has_analog[j] <- TRUE
          
          # Check if analog is within range
          future_lat <- point + displacements[j]
          within_range[j] <- (future_lat >= segments_display$start[i] && 
                                future_lat <= segments_display$end[i])
        }
      }
    }
    
    # Store center point displacement (for backward compatibility)
    center_idx <- which.min(abs(points - segments_display$center_lat[i]))
    if (length(center_idx) > 0) {
      segments_display$required_displacement[i] <- abs(displacements[center_idx])
    } else {
      segments_display$required_displacement[i] <- NA
    }
    
    # Count points with analogs
    segments_display$points_with_analog_anywhere[i] <- sum(has_analog)
    segments_display$points_with_analog_in_range[i] <- sum(within_range)
    
    # Calculate percentage of points with no analog within range
    if (segments_display$points_sampled[i] > 0) {
      segments_display$no_analog_within_range_pct[i] <- 100 - 
        (segments_display$points_with_analog_in_range[i] / 
           segments_display$points_sampled[i] * 100)
    } else {
      segments_display$no_analog_within_range_pct[i] <- 100  # Default to 100% if no points
    }
    
    # Find the minimum displacement among points that have analogs
    valid_displacements <- abs(displacements[has_analog])
    if (length(valid_displacements) > 0) {
      segments_display$min_required_displacement[i] <- min(valid_displacements, na.rm=TRUE)
      segments_display$no_analog[i] <- FALSE
    } else {
      segments_display$min_required_displacement[i] <- NA
      segments_display$no_analog[i] <- TRUE
    }
  }
  
  # Two criteria for complete habitat loss:
  # 1. Traditional: Minimum displacement exceeds segment size
  segments_display$complete_habitat_loss_traditional <- 
    segments_display$size <= segments_display$min_required_displacement
  
  # 2. New approach: 100% of points have no analog within range
  segments_display$complete_habitat_loss_new <- 
    segments_display$no_analog_within_range_pct >= 100
  
  # Combined approach (they should generally agree)
  segments_display$complete_habitat_loss <- 
    segments_display$complete_habitat_loss_traditional | 
    segments_display$complete_habitat_loss_new
  
  # Areas with no analog anywhere should be marked as complete habitat loss
  segments_display$complete_habitat_loss[segments_display$no_analog] <- TRUE
  
  # Existing display calculations
  segments_display$display_start <- ((segments_display$start - 0) / lat_range) * 100 - 50
  segments_display$display_end <- ((segments_display$end - 0) / lat_range) * 100 - 50
  segments_display$display_size <- segments_display$display_end - segments_display$display_start
  segments_display$display_center_lat <- ((segments_display$center_lat - 0) / lat_range) * 100 - 50
  
  scale_factor <- 100 / lat_range
  segments_display$display_displacement <- segments_display$required_displacement * scale_factor
  segments_display$display_min_displacement <- segments_display$min_required_displacement * scale_factor
  
  return(segments_display)
} 

calculate_min_range_size <- function(today_temp, future_temp, tolerance = 0.01) {
  # Find middle latitude (equator)
  min_lat <- min(today_temp$latitude)
  max_lat <- max(today_temp$latitude)
  middle_lat <- (min_lat + max_lat) / 2
  
  # Create results dataframe
  results <- data.frame(
    latitude = today_temp$latitude,
    temp_today = today_temp$temperature,
    distance_to_analog = NA,
    min_range_size = NA
  )
  
  # For each latitude in today's climate
  for (i in 1:nrow(results)) {
    current_lat <- results$latitude[i]
    current_temp <- results$temp_today[i]
    
    # Find all future latitudes with temperatures close to current temperature
    temp_diff <- abs(future_temp$temperature - current_temp)
    min_diff <- min(temp_diff)
    future_analog_indices <- which(abs(temp_diff - min_diff) < tolerance)
    future_analog_lats <- future_temp$latitude[future_analog_indices]
    
    # Determine which hemisphere we're in
    is_north_of_equator <- current_lat > middle_lat
    
    # Choose appropriate analog based on hemisphere
    if (is_north_of_equator) {
      # North of equator - prefer analogs north of equator but further south than current
      valid_analogs <- future_analog_lats[future_analog_lats <= current_lat & 
                                            future_analog_lats >= middle_lat]
      # If no valid analogs found, look for analogs further north
      if (length(valid_analogs) == 0) {
        valid_analogs <- future_analog_lats[future_analog_lats > current_lat]
      }
      # If still no valid analogs, consider all analogs
      if (length(valid_analogs) == 0) {
        valid_analogs <- future_analog_lats
      }
    } else {
      # South of equator - prefer analogs south of equator but further north than current
      valid_analogs <- future_analog_lats[future_analog_lats >= current_lat & 
                                            future_analog_lats <= middle_lat]
      # If no valid analogs found, look for analogs further south
      if (length(valid_analogs) == 0) {
        valid_analogs <- future_analog_lats[future_analog_lats < current_lat]
      }
      # If still no valid analogs, consider all analogs
      if (length(valid_analogs) == 0) {
        valid_analogs <- future_analog_lats
      }
    }
    
    # Calculate distances to future analogs
    distances <- abs(valid_analogs - current_lat)
    min_distance <- min(distances)
    results$distance_to_analog[i] <- min_distance
    
    # Calculate relative position on gradient (0 at equator, -1 at min lat, +1 at max lat)
    rel_position <- (current_lat - middle_lat) / (max_lat - middle_lat) * 2
    
    # Apply symmetrical scaling factor based on distance from equator
    # Use absolute value of rel_position to ensure symmetry
    factor <- 2 - abs(rel_position)
    factor <- max(factor, 1.0)
    
    # Calculate minimum range size with adjusted factor
    min_range_size <- factor * min_distance
    results$min_range_size[i] <- min_range_size
  }
  
  return(results)
}


## run

# Calculate latitudinal displacement
displacement <- calculate_latitudinal_displacement(today_temp, future_temp)

# Apply this to your segments
analyzed_segments <- find_critical_range_size(displacement, segments)

# Calculate the percentage of segments that will face complete habitat loss
percent_loss <- sum(analyzed_segments$complete_habitat_loss) / nrow(analyzed_segments) * 100
analyzed_segments$pct_domain <- analyzed_segments$size / 1000

combined_segment_results <- merge(results_sim, analyzed_segments, by = "segment_id", all.x = TRUE)

segments_lost <- subset(analyzed_segments, complete_habitat_loss == TRUE)

min_range <- calculate_min_range_size(today_temp1, future_temp1)
min_range$display_latitude <- -50 + 100 * (min_range$latitude  - 0) / 1000
min_range$display_range <- (min_range$min_range_size - 0) / 1000 * 100

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
max_critical_range <- max(min_range$display_range, na.rm = TRUE)

# plot range size vs habitat loss
plot3 <- ggplot(results_sim, aes(x = pct_domain*100, y = exceed_max_pct)) +
  annotate("rect", 
              xmin = min_critical_range, 
              xmax = max_critical_range, 
              ymin = -Inf, 
              ymax = Inf, 
              fill = "red", 
              alpha = 0.1) +
  annotate("rect", 
           xmin = 0, 
           xmax = min_critical_range, 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "red", 
           alpha = 0.3) +
  geom_pointdensity(size=0.8)+
  scale_color_gradient(low='#cbecff',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  # Add the cumulative species line
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
  labs(x = "Range size as\nproportion of domain", 
       y = "Area exceeding envelope (%)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  theme(legend.position="none",
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.line.y.right = element_line(color="blue"))

plot3_nored <- ggplot(results_sim, aes(x = pct_domain*100, y = exceed_max_pct)) +
  geom_pointdensity(size=0.8)+
  scale_color_gradient(low='#cbecff',high="#2171b5",trans="log10") + #scale_color for point, scale_fill for hex
  geom_smooth(method = "loess", span = 0.3,color = "black",se=F) +
  # Add the cumulative species line
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
  labs(x = "Latitudinal range \nsize (degrees)", 
       y = "Area exceeding envelope (%)")+
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20),
                     expand=c(0,0))+
  theme(legend.position="none",
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.line.y.right = element_line(color="blue"))

# Plot showing habitat loss by latitude
# Calculate percentage of segments at risk

plot_by_latitude <- ggplot(analyzed_segments, aes(x = display_center_lat, y = display_size)) +
  annotate("rect", 
           xmin = -Inf, 
           xmax = Inf, 
           ymin = min_critical_range, 
           ymax = max_critical_range, 
           fill = "red", 
           alpha = 0.1) +
  annotate("rect", 
           xmin = -Inf, 
           xmax = Inf, 
           ymin = 0, 
           ymax = min_critical_range, 
           fill = "red", 
           alpha = 0.3) +
  # Regular points
  geom_point(aes(color = complete_habitat_loss), size = 0.5,alpha = 0.3) +
  # No analog points 
  geom_point(data = subset(analyzed_segments, no_analog), 
             aes(x = display_center_lat, y = display_size),
             color = "red", size = 0.5, alpha = 0.3) +
  # distance to analog
  geom_line(data = displacement, aes(x = display_latitude, y = abs(display_displacement)), 
            color = "black", size = 1.5) +
  #minimum range size
  geom_line(data = min_range, 
            aes(x = display_latitude, y = display_range),
            color = "black", linetype = "dashed", size = 1.2) +
  scale_color_manual(values = c("blue", "red"), name = "Complete \nhabitat loss",
                     guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0),breaks=(seq(0, 100, 20)))+
  theme_classic() +
  labs(x = "Latitude", 
       y = "Range size \n(degrees latitude)") +
  theme(
    #legend.position="bottom",
    legend.position = c(0.98, 0.98),  # x and y coordinates (0-1 scale)
    legend.justification = c(1, 1),   # Anchor point at top right of legend box
    legend.background = element_rect(fill = "white", color = NA, linewidth = 0.5, linetype = "solid"),
    legend.margin = margin(3, 3, 3, 3),
    legend.spacing.y = unit(0.1, "cm"),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 9), 
    legend.text = element_text(size = 9))


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






fig2 <- plot_grid(plot1, plot4,plot2,plot3,plot_by_latitude,displace_loss,
          labels = c("A","B","C","D","E","F"), ncol = 3,nrow=2, align="hv")
ggsave(paste0(figures_folder, "fig2.png"), fig2, 
       width = 11, height = 6, dpi = 350)





















