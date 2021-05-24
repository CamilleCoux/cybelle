## extract and intgrate to main dataframe "track" environmental variables:
# - bathymetry from .asc file
# - chla from .nc file
# - sst from .nc file
# - site.to.coast (distance to nearest coast) calculated here
# output : m_grid_may20, an sf obect with the environmental variables averages over May 2020,
# plus the grid geometry and ID.

# load grid : 9852 cells that cover the whole Mediterranean basin
m_grid2 <- st_read("cybelle/data/med_grid.shp", crs=4326)

# to avoid confusions, change name of the grid to specify that it will contain the
# environmental measures for a given time. In this case, I chose mean values of over May 2020.
m_grid_may20 <- m_grid2

# rename variable to be consistent with the rest of the analysis
m_grid_may20 <- rename(m_grid_may20, site.to.coast = dtocoast)
# extract coordiantes from cell centroids, which are the points at which the environmental 
# variables will be extracted from 
centroids <- m_grid_may20 %>% 
  st_centroid %>%
  st_coordinates

# view grid
# mapview(centroids %>% 
#           as.data.frame %>% 
#           st_as_sf(coords=c("X", "Y"), crs=4326),
#         viewer.suppress = TRUE )


# bathymetry

bathymetry <- raster("cybelle/data/Bathymetry.asc")
m_grid_may20$bathymetry <- raster::extract(bathymetry, centroids)

# remove cells where bathymetry = NA, such that the grid is centered on the NW medtiterranean
# basin

m_grid_may20 <- m_grid_may20[-which(is.na(m_grid_may20$bathymetry)),]
mapview(m_grid_may20$geometry,
        viewer.suppress = TRUE )

# update centroids
centroids <- m_grid_may20 %>% 
  st_centroid %>%
  st_coordinates

# Monthly chla and sst data in .nc files

chla <- raster("cybelle/data/chla_july2020.nc")
m_grid_may20$chla <- raster::extract(chla, centroids)

sst <- raster("cybelle/data/sst_july2020.nc")
m_grid_may20$sst <- raster::extract(sst, centroids)


# Distance to nearest coast

# convert to sf object
centroids_sf <- 
  centroids %>% 
  as.data.frame %>% 
  st_as_sf(coords = c("X","Y"), crs=4326) 

# create the Mediterranean borders

# 1. 
med <- maps::map("world", 
                 fill=T, 
                 interior = T, 
                 xlim = c(-6, + 17), # select coordinates large enough to encompass
                 ylim = c(34, 47),   # the whole Mediterranean Sea
                 plot = F) 

# 2. turn into sf object with lines instead of points
med_sf <- maptools::map2SpatialPolygons(
  med, 
  ID=med$names, 
  proj4string = CRS("+init=epsg:4326")
  ) %>%
  st_as_sf %>%
  st_crop(xmin=-6, # crop to NW basin
          xmax=17, 
          ymin=34,
          ymax=47)

# 3. calculate distances from each observation data poit to the coast 
# (careful, with more cells this rapidly takes much longer to run)

site.to.coast <- med_sf %>% 
  st_distance(centroids_sf, by_element = F)

# 4. extract the minimal distance (i.e. nearest coast)
# head(site.to.coast) # don't do this

site.to.coast2 <- matrix(as.numeric(site.to.coast),
                   nrow = nrow(site.to.coast),
                   ncol = ncol(site.to.coast)) 
site.to.coastnum <- apply(site.to.coast2,2,min)/1000 # convert from m to km
head(site.to.coastnum)
range(site.to.coastnum)

m_grid_may20$site.to.coast <- site.to.coastnum

# clean workspace
rm( centroids, chla, site.to.coast, site.to.coast2, site.to.coastnum, m_grid2, med, med_sf, sst)

