#projection analysis - threshold exceedance and possible analogs

source("functions.R")

#open species envelope data
species_coord <- read.csv(paste0(output_directory,'envelopes/species_pixels_1982-1992.csv'))
species_coord$x <- ifelse(species_coord$x < 0, species_coord$x + 360, species_coord$x)
species_envelope <- read.csv(paste0(output_directory,'envelopes/climate_envelopes_1982-1992.csv'))



#run for both summers and annual means
#summer analysis tells us when coral will be bleaching all summer
#annual means tell us when coral will be under constant state of bleaching stress

mask <- rast(paste0(cmip_folder,'landseamask360l.tif'))
#fix mask - all latitudes above 84.5 are NA. make them 10
lat_values <- yFromCell(mask, 1:ncell(mask))
lat_raster <- setValues(mask, lat_values)
mask[lat_raster > 83.75 & is.na(mask)] <- 10

#open projected climate data
tos_summerN <- rast(paste0(cmip_folder,'tos/tos_OdecsummerN_modelmean_ssp585_201501-210012.tif'))
tos_summerS <- rast(paste0(cmip_folder,'tos/tos_OdecsummerS_modelmean_ssp585_201501-210012.tif'))
tos_lat_mask <- init(tos_summerN, fun= "y")
tos_summer <- ifel(tos_lat_mask >= 0, tos_summerN, tos_summerS)

ph_summerN <- rast(paste0(cmip_folder,'ph/ph_OdecsummerN_modelmean_ssp585_201501-210012.tif'))
ph_summerS <- rast(paste0(cmip_folder,'ph/ph_OdecsummerS_modelmean_ssp585_201501-210012.tif'))
ph_lat_mask <- init(ph_summerN, fun= "y")
ph_summer <- ifel(ph_lat_mask >= 0, ph_summerN, ph_summerS)


# harmonize coverage of climate data: exclude pixels where one variable is NA
global(tos_summer, "notNA")[[1]] #43976
global(tos_summer, "isNA")[[1]] #20824
global(ph_summer, "notNA")[[1]] #43597 
global(ph_summer, "isNA")[[1]] #21203
tos_summer[is.na(ph_summer)] <- NA
ph_summer[is.na(tos_summer)] <- NA
#after: 43158 notNA; 21642 isNA


scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
years <- c(2020,2030,2040,2050,2060,2070,2080,2090,2100)

#run analysis
for(i in 1:9){
  run_analog_analysis(
    species_list = species_coord,
    species_envelopes = species_envelope,
    future_temp_list = list(tos_summer[[i]]), 
    future_ph_list = list(ph_summer[[i]]),
    time_periods = c(years[[i]]),
    base_period = "1985",
    mask = mask,
    output_dir = output_directory,
    scenario = scenarios[[4]],
    model = "modelmean"
  )
  
  notifyitsdone(
    api_token = APItoken,
    channel = "coralspecies",
    event = "Analog Analysis Done",
    description = paste0("904 species analogs calculated for year ",years[[i]]),
    icon = "🔬"
  )
}



analogdual <- read.csv(paste0(output_directory,"dualanalog_modelmean_ssp585_1985-2100.csv"))
analogtos <- read.csv(paste0(output_directory, 'analog_tos_modelmean_ssp585_1985-2100.csv'))
analogph <- read.csv(paste0(output_directory, 'analog_ph_modelmean_ssp585_1985-2100.csv'))
analogdiff <- read.csv(paste0(output_directory,"analogdiff_modelmean_ssp585_1985-2100.csv"))

species_id <- 219556732
sp_coords <- species_coord[species_coord$species_id == species_id,]
sp_envelope <- species_envelope[species_envelope$id_no == species_id,]

analogdualhi <- analogdual[analogdual$species_id == species_id,]
analogtoshi <- analogtos[analogtos$species_id == species_id,]
analogphhi <- analogph[analogph$species_id == species_id,]
analogdiffhi <- analogdiff[analogdiff$species_id == species_id,]








### year of threshold exceedance per species
# for each species
# for each year
# calculate how many cells (and sum up the area) exceed envelope
# pvalf_temp > pmax_temp. pvalf_ph < pmin_ph

# record difference (how much future condition exceeds envelopes)
# pvalf_temp - pmax_temp. pvalf_ph - pmin_ph
# record year when all cells in distribution exceed envelope for temp and ph

# possible analogs
# mean distance to analog per species
# for dual, univariate, and diff
# are any analogs within current distribution



# are there possible analogs? where? when do analogs disappear? first loss of suitablility and last suitable period
