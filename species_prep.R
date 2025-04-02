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





