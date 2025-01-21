# figures
library(cowplot)

source("functions.R")

#open summary data
coralsp126 <- fread(file=paste0(output_directory,"coralsp_summary_",scenarios[[1]],"_1985-2100.csv"))
coralsp245 <- fread(file=paste0(output_directory,"coralsp_summary_",scenarios[[2]],"_1985-2100.csv"))
coralsp370 <- fread(file=paste0(output_directory,"coralsp_summary_",scenarios[[3]],"_1985-2100.csv"))
coralsp585 <- fread(file=paste0(output_directory,"coralsp_summary_",scenarios[[4]],"_1985-2100.csv"))


# Winners and losers - table showing which species have later/sooner threshold exceedance dates
# Current range area - indicator of if it’s an endemic species?






# Number of species that have exceeded thresholds under each RCP





#percentage of area per species that remains suitable refugia over time and environmental conditions at those sites
p1 <- ggplot(coralsp585) +
  geom_line(aes(x = year, 
                 y = pctrefugia_temp,
                group = species_id, 
                 color = pvalfExceedph_temp), # how much pH differs at temp refugia
             alpha = 0.6) +
  scale_color_gradient2(midpoint = 0,
                        low = "orange", 
                        mid = "white") +
  labs(x = "Year",
       y = "Percent area of temperature refugia",
       color = "pH difference at refugia") +
  theme_minimal()

# Plot 2: pH refugia and their temperature conditions  
p2 <- ggplot(coralsp585) +
  geom_line(aes(x = year, 
                 y = pctrefugia_ph,
                 group = species_id, 
                 color = pvalfExceedtemp_ph), # how much temp differs at pH refugia
             alpha = 0.6) +
  scale_color_gradient2(midpoint = 0,
                        low = "blue", 
                        mid = "white",
                        high = "red") +
  labs(x = "Year",
       y = "Percent area of pH refugia",
       color = "Temperature difference at refugia") +
  theme_minimal()
plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)







#difference in conditions at current location
#Quantification of how different future climate will be - a measure of “adaptation potential” or “resilience requirement”


# Historic vs current vs future mean temp. Historic vs current vs future summer temp

# Distance from current range to closest dual and individual analogs



# rejoin with original shapefile - where are hotspots of of multiple species exceding thresholds? 

































# other plots

# scatter plot year when each species complete distribution area exceeds envelopes - for temperature, ph, and both
ggplot(coralsp585) +
  geom_point(aes(x = year, y = totalarea, color = "Temperature"), 
             data = . %>% filter(!is.na(completeExceed_temp)), alpha = 0.6) +
  geom_point(aes(x = year, y = totalarea, color = "pH"),
             data = . %>% filter(!is.na(completeExceed_ph)), alpha = 0.6) +
  geom_point(aes(x = year, y = totalarea, color = "Both"),
             data = . %>% filter(!is.na(completeExceed_tempph)), alpha = 0.6) +
  scale_color_manual(values = c("Temperature" = "blue", 
                                "pH" = "purple", 
                                "Both" = "green")) +
  labs(title = "Complete Loss of Suitable Habitat Over Time",
       x = "Year",
       y = "Total Area",
       color = "Loss Type") +
  theme_minimal()
