library(sf)
library(terra)

source("config_dir.R")

# 904 species of reef-building coral polygons split in 3 shapefiles
coral1 <- vect(paste0(iucn_sp_folder,"REEF_FORMING_CORALS_PART1.shp"))
coral2 <- vect(paste0(iucn_sp_folder,"REEF_FORMING_CORALS_PART2.shp"))
coral3 <- vect(paste0(iucn_sp_folder,"REEF_FORMING_CORALS_PART3.shp"))

coral <- rbind(coral1, coral2, coral3)

#shape_area needs to be recalculated according to metadata
coral$distribution_area_km2 <- expanse(coral, unit="km")

writeVector(coral, paste0(iucn_sp_folder,"coral_all.gpkg"),driver="GPKG", overwrite=TRUE) 
rm(coral)  
gc()  








#focal stats for increased coastal coverage
scenarios <- c("ssp126","ssp245","ssp370","ssp585") 
variables <- c("tos","ph")
models <- c('modelmean','modelmedian','modelmeanfoc')

sst <- rast(paste0(cmip_folder,'tos/tos_Omon_modelmean_hist_185001-201412.tif'))
focal(sst, w=matrix(1, 3, 3),fun=mean,na.rm=T,NAonly=T,
      filename=paste0(cmip_folder,'tos/tos_Omon_modelmeanfoc_hist_185001-201412.tif'))
ph <- rast(paste0(cmip_folder,'ph/ph_Omon_modelmean_hist_185001-201412.tif'))
focal(ph, w=matrix(1, 3, 3),fun=mean,na.rm=T,NAonly=T,
      filename=paste0(cmip_folder,'ph/ph_Omon_modelmeanfoc_hist_185001-201412.tif'))

for(scenario in scenarios){
  for(var in variables){
    vari <- rast(paste0(cmip_folder,var,'/',var,'_Omon_',models[1],'_',scenario,'_201501-210012.tif'))
    focal(vari, w=matrix(1, 3, 3),fun=mean,na.rm=T,NAonly=T,
          filename=paste0(cmip_folder,var,'/',var,'_Omon_',models[3],'_',scenario,'_201501-210012.tif'))
  }
}




#add on omega arag

#data cleaning: imports on a 0.5 360.5 0.5 180.5 grid, not aligning with the rest of the data
#in cdo: cdo setgrid,r360x180 Aragonite_median_historical.nc Aragonite_median_historical_fix.nc
#need to shift 20 degrees east - visibly misaligned with the rest of datasets
custom_rotate <- function(rast_obj, pivot_lon = 0) {
  # Get original extent
  e <- ext(rast_obj)
  
  # Determine which parts to split
  right_extent <- c(pivot_lon, e$xmax, e$ymin, e$ymax)
  left_extent <- c(e$xmin, pivot_lon, e$ymin, e$ymax)
  
  # Crop the two parts
  right_part <- crop(rast_obj, ext(right_extent))
  left_part <- crop(rast_obj, ext(left_extent))
  
  # Calculate shift amounts
  shift_right <- e$xmin - pivot_lon
  shift_left <- e$xmax - pivot_lon
  
  # Apply shifts
  right_shifted <- terra::shift(right_part, dx = shift_right)
  left_shifted <- terra::shift(left_part, dx = shift_left)
  
  # Merge the parts
  result <- merge(right_shifted, left_shifted)
  
  # Ensure correct extent (might need adjustment based on your specific case)
  ext(result) <- c(e$xmin, e$xmax, e$ymin, e$ymax)
  
  return(result)
}

arag_hist <- rast(paste0(cmip_folder, "arag/Aragonite_median_historical_fix.nc"))
arag_shifted <- custom_rotate(arag_hist,pivot_lon=339)
arag_shifted <- ifel(arag_shifted == 0, NA, arag_shifted)

arag_hist_foc <- focal(arag_shifted, w=matrix(1, 3, 3), fun=mean, na.rm=TRUE, NAonly=TRUE)
writeRaster(arag_hist_foc, 
            filename=paste0(cmip_folder, "arag/Aragonite_median_historical_foc.nc"),
            overwrite=TRUE)


for(scen_file in c("126", "245", "370", "585")) {
  # Read in the file
  arag_scen <- rast(paste0(cmip_folder, "arag/Aragonite_median_ssp", scen_file, "_fix.nc"))
  
  arag_shifted <- custom_rotate(arag_scen,pivot_lon=339)
  arag_shifted <- ifel(arag_shifted == 0, NA, arag_shifted)
  
  # Apply focal stats
  focal(arag_shifted, w=matrix(1, 3, 3), fun=mean, na.rm=T, NAonly=T,
        filename=paste0(cmip_folder, "arag/Aragonite_median_ssp", scen_file, "_foc.nc"),overwrite=T)
}







plot(sst[[1]])
lines(world$long, world$lat, col = "black")

plot(arag_shifted[[1]])
lines(world$long, world$lat, col = "black")
