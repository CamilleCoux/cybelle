# Code by Camille Coux
# 


library(magrittr)
library(tidyr)
library(dplyr)
library(unmarked)
library(ggplot2)
library(sf)
library(raster)
library(ncdf4)
library(mapview)
library(maps)
library(maptools)


options(scipen = 999)
theme_set(theme_minimal())

# download data (ask me by writing at camille.coux@orange.fr)

# import obs data
track <- read.csv2("data/track.csv", row.names = 1)

# import grid:
m_grid <- st_read("data/med_grid.shp", crs=4326)

# view grid in browser:
# mapview(m_grid$geometry,viewer.suppress = TRUE )

# intersect obs data with grid cells to get the grid_id column in the track dataset :
inter <- track %>%
  st_as_sf(coords=c("long", "lat"), crs=4326) %>%
  st_intersection(m_grid, track_sf) %>%
  rename(grid_id = FID)
track <- left_join(track, inter%>%dplyr::select(grid_id, index))

# select observations for a given species, e.g. striped dolphin species, 
# i.e. "dauphin bleu et blanc" in french:
presences <- track[grep("Dauphin Bleu", track$species),]
# selection of all absences  
absences <- track[which(is.na(track$species)),]

# merge
striped_dolphin <- rbind(presences, absences)

# takeout presence only obseravtions
striped_dolphin <- striped_dolphin[-grep("ponctuelle", striped_dolphin$protocole),]

# keep observations from the summer months only since they concentrate most obseravtions
presences$month %>% table %>% barplot 
striped_dolphin <- striped_dolphin %>%
  filter(month %in% 6:9)

# check which years have enough data : keep data from 2015 to 2020
table(striped_dolphin$month, striped_dolphin$year)
striped_dolphin <- striped_dolphin %>% filter(year %in% c(2015:2020))

# select columns necessary for the analysis, and remove NAs
d <- striped_dolphin %>%
  dplyr::select(bathymetry, site.to.coast, chla, sst) %>%
  complete.cases
striped_dolphin <- striped_dolphin[d,]

# format data to create the unmarked data frame object
# 1. create An RxJ matrix of the detection, non-detection data, where: 
# R =  number of sites
# J = maximum number of sampling periods per site

RJ_mat <- striped_dolphin %>% 
  group_by(grid_id, year) %>%
   summarise(n=sum(n, na.rm = T)) %>%
   pivot_wider(names_from = year, values_from=n, names_prefix = "Y")
RJ_mat <- RJ_mat[,c(1, 7, 3, 2, 6, 4, 5)]

# Select observation variables : chla and SST
# calculate the mean chla and sst values for each grid cell:
chla <- striped_dolphin %>% 
  group_by(grid_id, year) %>%
  summarise(chla = mean(chla, na.rm = T)) %>%
  pivot_wider(names_from = year, values_from=chla, names_prefix = "chla") %>%
  as.data.frame
chla <- chla[,c(1, 7, 3, 2, 6, 4, 5)]


sst <- striped_dolphin %>% 
  group_by(grid_id, year) %>%
  summarise(sst = mean(sst, na.rm = T)) %>%
  pivot_wider(names_from = year, values_from=sst, names_prefix = "sst") %>%
  as.data.frame
sst <- sst[,c(1, 7, 3, 2, 6, 4, 5)]

# merge chla and sst into single observation matrix:
obscovs <- list(chla[,-1] , sst[,-1])
names(obscovs) <- c("chla", "sst")

# Select site covariates : bathymetry and distance from site to coast
# calculate the mean bathymetry and site.to.coast values for each grid cell:
sitecovs <- striped_dolphin %>% 
  group_by(grid_id) %>%
  summarise(bathymetry.sites = round(mean(bathymetry, na.rm=T), digits = 1),
            site.to.coast = round(mean(site.to.coast, na.rm=T), digits = 1)) 

umf <- unmarkedFrameOccu(y = RJ_mat[,-1] %>% as.data.frame, 
                         siteCovs = sitecovs[,-1] %>% as.data.frame, 
                         obsCovs = obscovs)


# scale covariates and store values for later
sc <- scale(siteCovs(umf))
siteCovs(umf) <- sc
scobs <- scale(obsCovs(umf))
obsCovs(umf) <- scobs
umf                     # look at data
summary(umf)            # summarize      

occu.model = occu(~chla + sst ~ bathymetry.sites+site.to.coast, umf)



#--------- carte sans model-averaging (debut) --------------------------
# read in Med grid with environmental variables extracted for May 2020:
source("code/extract_env_data.R")
# there will be warnings, they're ok for this purpose

# check NAs
m_grid_may20 %>% apply(., 2, is.na) %>% colSums

# remove lines with NAs:
m_grid_may20 <- m_grid_may20[-which(is.na(m_grid_may20$chla)), ]

# values used to standardise the unmarked dataframe variables: need to apply
# the same standardisation values to all cells of m_grid_may20 
sc
mean_Bathy <- attributes(sc)$`scaled:center`[1]
sd_Bathy <- attributes(sc)$`scaled:scale`[1]
bathy.s <- (m_grid_may20$bathymetry - mean_Bathy) / sd_Bathy

# occupancy model estimates for bathymetry
summary(occu.model)
(beta <- coef(occu.model, type="state"))
logit.psi <- beta[1] + beta[2]*bathy.s 
psi <- exp(logit.psi) / (1 + exp(logit.psi))

# get quantiles of occupancy estimates
grid_occ <- quantile(psi,probs=seq(0, 1, 0.1), na.rm=T) 
round(grid_occ,2)

# maps

# 1. the actual bathymetry measures
mapview(m_grid_may20, viewer.suppress=T, zcol="bathymetry")

# 2. occupancy predictions based on bathymetry variable
labs <- c("0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "90-100%")
m_grid_may20$bathy_prediction <- psi
mapview(m_grid_may20, viewer.suppress=T, zcol="bathy_prediction")

# 3. occupancy predictions based on bathymetry variable - discretised using quantiles
m_grid_may20$psi_bathy_quantiles <- cut(psi,breaks= grid_occ,labels=labs) 
mapview(m_grid_may20, viewer.suppress=T, zcol="psi_bathy_quantiles")


# The distance to the coast wasn't significant

# What about the chla effect ?
scobs
mean_chla =attributes(scobs)$`scaled:center`[1]
sd_chla= attributes(scobs)$`scaled:scale`[1]
chla.s <- (m_grid_may20$chla - mean_chla) / sd_chla


# occupancy estimates from model p(chla)
(beta.det <- coef(occu.model, type="det"))
logit.p <- beta.det[1] + beta.det[2]*chla.s 
p <- exp(logit.p) / (1 + exp(logit.p))


# get quantiles of occupancy estimates
grid_occ_p <- quantile(p,probs=seq(0, 1, 0.1), na.rm=T) 
round(grid_occ_p,2)

# maps

# 1. the actual chla measures
mapview(m_grid_may20, viewer.suppress=T, zcol="chla")

# 2. actual chla measures but binned in quantiles
grid_chla_quantiles <- quantile(m_grid_may20$chla, probs=seq(0, 1, 0.1), na.rm=T)
m_grid_may20$chla_quantiles <- cut(m_grid_may20$chla,breaks= grid_chla_quantiles,labels=labs) 
mapview(m_grid_may20, viewer.suppress=T, zcol="chla_quantiles")

# 3. occupancy predictions based on chla variable
m_grid_may20$chla_prediction <- p
fig <- mapview(m_grid_may20, viewer.suppress=T, zcol="chla_prediction")

# To make these into .png files, this used to work, but generates an error now.
# Maybe because of my recent R update ?
# mapshot(m, file ="striped_dolphin_prediction_chla.png", 
# map.types = "CartoDB.Positron")


# 4. occupancy predictions based on chla variable - discretised using quantiles
m_grid_may20$p_chla_quantiles <- cut(p,breaks= grid_occ_p,labels=labs) 
mapview(m_grid_may20, viewer.suppress=T, zcol="p_chla_quantiles")


