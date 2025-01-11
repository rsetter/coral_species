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


