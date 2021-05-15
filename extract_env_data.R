## import et export des fichiers .nc

options(scipen = 999)
# crsproj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

library(magrittr)
library(raster)
library(lubridate)
library(dplyr)
library(sf)
library(rgdal)
library(ncdf4)
library(geosphere)
library(ggplot2)
# Import OBS DATA
track <- read.csv2("data/tout par annees/avec_trackings_mails/track_2020-niveau-proto-abond.csv", row.names = 1)

# load grid data
m_grid <- st_read("data/env_data/mediterranee shapefiles/toutemed_grid.shp", crs=4326)
centroids <- st_transform(m_grid, 4326) %>% 
  st_centroid %>%
  st_coordinates

# bathymetrie

bathymetrie <- raster("data/env_data/bathymetrie/Bathy_CC.asc")
m_grid$bathymetry <- raster::extract(bathymetrie, centroids)

#################################################################################


# Monthly chla data in .nc files, that are listed by year and julian years
# I need the mean for July 2020.
# To find out which file to open and figure out which are the corresponding julian days:
tmp <- as.POSIXlt("31Jul20", format = "%d%b%y")
tmp$yday

files=list.files(paste("/Users/camillecoux/Documents/Cybelle Planete/CybelleMed/data/env_data/chla/monthly/A20201832020213.L3m_MO_CHL.x_chlor_a.nc")) 
obsXY <- m_grid %>% dplyr::select(long, lat)

track$chla <- raster::extract(chlaj, obsXY)



j=2020
files=list.files(paste("/Users/camillecoux/Documents/Cybelle Planete/CybelleMed/data/env_data/chla/weekly/chla_8D_2020")) %>% 
  grep(j, ., value = T)
# stack up files per year
chlaj <- raster::stack(paste("/Users/camillecoux/Documents/Cybelle Planete/CybelleMed/data/env_data/chla/weekly/chla_8D_2020/", files, sep="")) 
# datapoints matrix:
obsXY <- track %>% dplyr::select(Evenement.Longitude, Evenement.Latitude)
track$chla <- raster::extract(chlaj, obsXY)

chla <- data.frame(grid.id = m_grid$grid_id, jun=NA, jul=NA, aug=NA, sep=NA)
for(m in 1:4){
  chla[,m+1] <- raster::extract(chlaj[[1]], centroids)
}






#################################################################################

# # pour la SST: data organisées par numéros de semaines depuis le début de chaque année.
# # je vais extraire la 1ère des data de tracking : c'est la semaine 26 donc la 26e de la
# # rasterbrick de 2009:
# j=2019
# files= list.files(paste("../env_data/sst/", j, sep="") ,pattern='*.nc$',full.names=TRUE)
# files <- files[grep("NSST", files)] # select the night temperatures
# sstj <- stack(files) # stack up files per year
# sstj <- crop(sstj, e) # crop the global data to just the Mediterranean
# names(sstj) %<>% gsub("[A-z]...", "", .) # simplify names to the 8 day bin number
# sstj[[26]]
# m_grid$sst <- raster::extract(sstj[[26]], centroids)
# rm(sstj)

j=2019
files=list.files(paste("../../../Dropbox/cybelle/env_data/sst/monthly")) %>% grep(j, ., value = T)
# pour les mois de juin à septembre: mois 6 à 9
files <- files[6:9]
# stack up files per year
sstj <- raster::stack(paste("../../../Dropbox/cybelle/env_data/sst/monthly/", files, sep="")) 
# make a df with sites as rows and monthly chla values as cols
sst <- data.frame(grid.id = m_grid$grid_id, jun=NA, jul=NA, aug=NA, sep=NA)
for(m in 1:4){
  sst[,m+1] <- raster::extract(sstj[[1]], centroids)
}


# Effort measure:  Calculate legth of trajectories in each grid cell for each year
##################################################################################
