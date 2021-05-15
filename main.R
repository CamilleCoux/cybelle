# script based on dauphin_occ2

library(magrittr)
library(tidyr)
library(dplyr)
library(AICcmodavg)
library(unmarked)
library(ggplot2)
library(sf)
library(mapview)
library(MuMIn)
library(RColorBrewer)


theme_set(theme_minimal())

# import obs data
track <- read.csv2("/Users/camillecoux/Dropbox/cybelle/data/track_2009_2020_simple.csv", row.names = 1)

# import grid:
m_grid <- st_read("data/env_data/mediterranee shapefiles/toutemed_grid.shp", crs=4326)
# mapview(m_grid$geometry,viewer.suppress = TRUE )


stripeddolphin_sf <- striped_dolphin %>%
  filter(n>0) %>%
  st_as_sf(coords=c("long", "lat"), crs=4326)
mapview(stripeddolphin_sf, zcol="protocole2", viewer.suppress = TRUE )

mapview(., zcol="protocole2", col=mycol[3], col.regions=mycol[3], cex=1  , legend=F) +
  
stripeddolphin_sf$geometry

# intersect obs data with grid cells:

inter <- track %>%
  st_as_sf(coords=c("long", "lat"), crs=4326) %>%
  st_intersection(m_grid, track_sf) %>%
  rename(grid_id = FID)
track <- left_join(track, inter%>%dplyr::select(grid_id, index))


# select observations for the striped dolphin species, i.e. "dauphin bleu et blanc" in french:
presences <- track[grep("Dauphin Bleu", track$species),]
# selection of all absences  
absences <- track[which(is.na(track$species)),]
# write.csv2(absences, "../cartes/absences.csv ")

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

d <- striped_dolphin %>%
  select(grid_id, year, n)


d <- striped_dolphin %>%
  select(bathymetry.sites, site.to.coast, chla, sst) %>%
  complete.cases
striped_dolphin <- striped_dolphin[d,]

# create An RxJ matrix of the detection, non-detection data, where: 
# R =  number of sites
# J = maximum number of sampling periods per site

RJ_mat <- striped_dolphin %>% 
  group_by(grid_id, year) %>%
   summarise(n=sum(n, na.rm = T)) %>%
   pivot_wider(names_from = year, values_from=n, names_prefix = "Y")
RJ_mat <- RJ_mat[,c(1, 7, 3, 2, 6, 4, 5)]

RJ_mat[,2:7] %>% is.na %>% rowSums %>% sort %>% unique


# Select observation variables : chla and SST
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

# merge chla and sst into sing observation matrix:
obscovs <- list(chla[,-1] , sst[,-1])
names(obscovs) <- c("chla", "sst")

# Select site covariates : bathymetry and distance from site to coast
# bathymetry: need to calculate means for each site based on the bathymetry
# values from each observation:
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

models <- dredge(occu.model)
#models with delta.aicc < 4
summary(model.avg(models, subset = delta < 2))




#--------- carte sans model-averaging (debut) --------------------------
# read in Med grid with bathy mean values
m_grid %>% apply(., 2, is.na) %>% colSums
# i need the mean bathys organised by grid_id. 
# juste pour tester un exemple j'ai pris les valeurs pour la toute première ligne de 2009
# pour la bathy mais aussi pour la chla et la sst

# valeurs utilis?es pour standardiser le tableau de covariables 
# (voir plus haut, avant l'ajustement des mod?les)
sc
mean_Bathy = 1437.83714
sd_Bathy = 1013.29677
bathy.s <- (m_grid$bathymetry - mean_Bathy) / sd_Bathy

(beta <- coef(occu.model, type="state"))
logit.psi <- beta[1] + beta[2]*bathy.s 
psi <- exp(logit.psi) / (1 + exp(logit.psi))


# occupancy estimates from model psi(bathy)
#exemple modele 20
summary(occu.model)

occ_matrix <- as.vector(psi) 

#m_grid@data$estim_occ <- occ_matrix
#writeOGR(obj= m_grid, dsn=".", layer="m_grid2", driver="ESRI Shapefile")

# get quantiles of occupancy estimates
grid_occ <- quantile(occ_matrix,probs=seq(0, 1, 0.11), na.rm=T) 
round(grid_occ,2)

# get quantiles of bathymetry
grid_bathy <- quantile((m_grid$bathymetry - mean_Bathy)/sd_Bathy,probs=seq(0, 1, 0.11), na.rm=T)
grid_bathy

# palette of greys (9 levels)
cols <- brewer.pal(n=9,name="Greys") 

# convert to spatialpolydataframe
m_grid2 <- as(m_grid, "Spatial")
# map bathy and occupancy estimates side by side
par(mfrow=c(1,2))

# Bathy cover first
lcols_f <- cut((m_grid2$bathymetry - mean_Bathy)/sd_Bathy,breaks= grid_bathy,labels=cols) # discretize bathy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='bathymétrie') # background grid grille de fond
plot(m_grid2, col=as.character(lcols_f), add=TRUE) # add bathy

# occupancy estimates second
lcols_psi <- cut(occ_matrix,breaks= grid_occ,labels=cols) # discretize occupancy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='occupancy') # background grid
plot(m_grid2, col=as.character(lcols_psi), add=TRUE) # add occupancy estimates

# donc là les deux plots sont identiques: l'occupancy est directement corrélée avec la bathy.


# maintenant je veux faire pareil mais avec les variables de détection: chla et effort.km
# valeurs utilis?es pour standardiser le tableau de covariables 
# (voir plus haut, avant l'ajustement des mod?les)
scobs
mean_chla = 0.2020719
sd_chla= 0.05781726
chla.s <- (m_grid$chla - mean_chla) / sd_chla

mean_eff <- 4.6059793
sd_eff <- 3.60565977
eff.s <- (m_grid$dtocost - mean_eff) / sd_eff


(beta.det <- coef(fm20, type="det"))
logit.p <- beta.det[1] + beta.det[2]*chla.s + beta.det[3]*eff.s
p <- exp(logit.p) / (1 + exp(logit.p))


# occupancy estimates from model p(chla)


occ_matrix <- as.vector(p) 

# get quantiles of occupancy estimates
grid_occ <- quantile(occ_matrix,probs=seq(0, 1, 0.11), na.rm=T) 
round(grid_occ,2)

# get quantiles of chla and effort
grid_chla <- quantile((m_grid$chla - mean_chla)/sd_chla,probs=seq(0, 1, 0.11), na.rm=T)
grid_chla
# palette of greys (9 levels)
cols <- brewer.pal(n=6,name="Greys") # 6 because there are 6 different values in grid.occ before = 1

# convert to spatialpolydataframe
m_grid2 <- as(m_grid, "Spatial")
# map bathy and occupancy estimates side by side
par(mfrow=c(1,2))

# chla cover first
lcols_f <- cut((m_grid2$chla - mean_chla)/sd_chla,breaks= grid_chla,labels=cols) # discretize bathy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='chla') # background grid grille de fond
plot(m_grid2, col=as.character(lcols_f), add=TRUE) # add bathy

# detection estimates second
lcols_p <- cut(occ_matrix,breaks= grid_occ[1:7],labels=cols) # discretize occupancy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='detection') # background grid
plot(m_grid2, col=as.character(lcols_p), add=TRUE) # add occupancy estimates


grid_effort <- quantile((m_grid$dtocost - mean_eff)/sd_eff, probs=seq(0, 1, 0.11), na.rm=T)


#--------- carte sans model-averaging (fin) --------------------------

#--------- carte avec model-averaging (début)-------------------------
# model selection
# je ne vais prendre que ceux qui ont un delat AIC inférieur à 4
as.matrix(cbind(mask4,labels[mask4])) 
tri <- cbind(labels[mask4],round(aic[mask4],2),round(aic[mask4]-aic[mask4[1]],2))
tri[which(as.numeric(tri[,3]) < 10),1]
fms <- fitList(m1=fm20,m2=fm21,m3=fm22,m4=fm26,m5=fm27,m6=fm28)
modSel(fms)

# model-averaging of psi for a grid of values for Bathy, and mean value for dtocoast
new_bathy <- (m_grid2@data$bathymetry - mean_Bathy)/sd_Bathy
# remove the NAs from bathy:
new_bathy <- new_bathy[-which(is.na(new_bathy))]
newdata <- data.frame(bathymetry=new_bathy,dtocoast=0)
# newdata <- newdata[-which(is.na(newdata$bathymetry)),]
# newdata <- as.matrix(newdata)
occ_pred_res <- predict(fms,type='state',newdata=newdata, appendData =TRUE, )
occ_pred <- occ_pred_res$Predicted
occ_fit <- 1/(1+exp(-occ_pred)) # back-transformed
#m_grid@data$estim_occ <- occ_matrix
#writeOGR(obj= m_grid, dsn=".", layer="m_grid2", driver="ESRI Shapefile")

# get quantiles of occupancy estimates
grid_occ <- quantile(occ_pred,probs=seq(0, 1, 0.11), na.rm=T) 
round(grid_occ,2)

# get quantiles of bathymetry
grid_bathy <- quantile((m_grid$bathymetry - mean_Bathy)/sd_Bathy,probs=seq(0, 1, 0.11), na.rm=T)
grid_bathy

# palette of greys (9 levels)
cols <- brewer.pal(n=9,name="Greys") 

# convert to spatialpolydataframe
m_grid2 <- as(m_grid, "Spatial")
# map bathy and occupancy estimates side by side
par(mfrow=c(1,2))

# Bathy cover first
lcols_f <- cut((m_grid2$bathymetry - mean_Bathy)/sd_Bathy,breaks= grid_bathy,labels=cols) # discretize bathy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='bathymétrie') # background grid grille de fond
plot(m_grid2, col=as.character(lcols_f), add=TRUE) # add bathy

# occupancy estimates second
lcols_psi <- cut(occ_pred,breaks= grid_occ,labels=cols) # discretize occupancy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='occupancy') # background grid
plot(m_grid2, col=as.character(lcols_psi), add=TRUE) # add occupancy estimates

# donc là les deux plots sont identiques: l'occupancy est directement corrélée avec la bathy.




#--------- carte avec model-averaging (fin)-------------------------

#--------- model averaging (debut) -----------------------------------

# model selection
# je ne vais prendre que ceux qui ont un delat AIC inférieur à 4
as.matrix(cbind(mask4,labels[mask4])) 
tri <- cbind(labels[mask4],round(aic[mask4],2),round(aic[mask4]-aic[mask4[1]],2))
tri[which(as.numeric(tri[,3]) < 10),1]
fms <- fitList(m1=fm20,m2=fm21,m3=fm22,m4=fm26,m5=fm27,m6=fm28)
modSel(fms)

# model-averaging of psi for a grid of values for Bathy, and mean value for SST and Chla
new_bathy <- (m_grid2@data$bathymetry - mean_Bathy)/sd_Bathy
# remove the NAs from bathy:
new_bathy <- new_bathy[-which(is.na(new_bathy))]
newdata <- data.frame(bathymetry=new_bathy,chla=0,sst=0, dtocoast=0)
# newdata <- newdata[-which(is.na(newdata$bathymetry)),]
# newdata <- as.matrix(newdata)
occ_pred_res <- predict(fms,type='state',newdata=newdata, appendData =TRUE)
occ_pred <- occ_pred_res$Predicted
# get quantiles of occupancy estimates
grid_occ <- quantile(occ_pred,probs=seq(0, 1, 0.11))
round(grid_occ,2)
# get quantiles of bathymetry
grid_bathy <- quantile((m_grid2@data$bathymetry - mean_Bathy)/sd_Bathy,probs=seq(0, 1, 0.11), na.rm=T)
grid_bathy

# palette of greys (9 levels)
cols <- brewer.pal(n=9,name="Greys") 

# map bathy and occupancy estimates side by side
par(mfrow=c(1,2))

# Bathy cover first
lcols_f <- cut((m_grid2$bathymetry - mean_Bathy)/sd_Bathy,breaks= grid_bathy,labels=cols) # discretize bathy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='bathymétrie') # background grid grille de fond
plot(m_grid2, col=as.character(lcols_f), add=TRUE) # add bathy

# occupancy estimates second
lcols_psi <- cut(occ_pred,breaks= grid_occ,labels=cols) # discretize occupancy using quantiles
plot(m_grid2, axes=TRUE,border="gray",main='model-averaged estimated occupancy') # background grid
plot(m_grid2, col=as.character(lcols_psi), add=TRUE) # add occupancy estimates

#--------- model averaging (fin) -----------------------------------


#--------- relation entre occupancy et bathy, sst et chla - model averaging (debut) -----------------------------------

# model selection
# je ne vais prendre que ceux qui ont un delat AIC inférieur à 4
as.matrix(cbind(mask4,labels[mask4])) 
tri <- cbind(labels[mask4],round(aic[mask4],2),round(aic[mask4]-aic[mask4[1]],2))
tri[which(as.numeric(tri[,3]) < 4),1]
fms <- fitList(m1=fm20,m2=fm21,m3=fm22,m4=fm26,m5=fm27,m6=fm28)
modSel(fms)


mean_cov <- apply(sitecovs,2,mean)
sd_cov <- apply(sitecovs,2,sd)
min_cov <- apply(sitecovs,2,min)
max_cov <- apply(sitecovs,2,max)

# model-averaging of psi for a grid of values for Bathy, 
# and mean value for dtocoast
grid_bathy_nonst <- seq(min_cov[1],max_cov[1],length=1000)
grid_bathy <- (grid_bathy_nonst - mean_cov[1])/sd_cov[1]
newdata <- data.frame(bathymetry=grid_bathy,dtocoast=0)
occ_pred_res <- predict(fms,type='state',newdata=newdata, appendData =TRUE)
occ_pred_bathy <- occ_pred_res$Predicted

# bathy effect on occupancy
plot(grid_bathy_nonst, occ_pred_bathy, 
     xlab= "Bathymétrie (m)", ylab="Probabilité d'occupation", 
     main="Relation entre occupation et bathymétrie", 
     ylim=c(0,1),las=1,type='n',lwd=3)
polygon(c(rev(grid_bathy_nonst), grid_bathy_nonst),
        c(rev(occ_pred_res$upper), occ_pred_res$lower), 
        col = 'grey80', border = NA)
lines(grid_bathy_nonst, occ_pred_bathy, col="blue",ylim=c(0,1),las=1,type='l',lwd=3)

# model-averaging of psi for a grid of values for dtocoast, 
# and mean value for bathy and Chla
grid_dtocoast_nonst <- seq(min_cov[2],max_cov[2],length=1000)
grid_dtocoast <- (grid_dtocoast_nonst - mean_cov[2])/sd_cov[2]
newdata <- data.frame(dtocoast=grid_dtocoast,bathymetry=0)
occ_pred_res <- predict(fms,type='state',newdata=newdata, appendData =TRUE)
occ_pred_dtocoast <- occ_pred_res$Predicted
# dtocoast effect on occupancy
plot(grid_dtocoast_nonst, occ_pred_dtocoast, xlab= "Dist à la côte", ylab="Probabilité d'occupation", 
     main="Relation entre occupation et \n distance à la côte", ylim=c(0,1),las=1,type='n',lwd=3)
polygon(c(rev(grid_dtocoast_nonst), grid_dtocoast_nonst), c(rev(occ_pred_res$upper), occ_pred_res$lower), col = 'grey80', border = NA)
lines(grid_dtocoast_nonst, occ_pred_dtocoast, col="blue",ylim=c(0,1),las=1,type='l',lwd=3)


#--------- relation entre occupancy et bathy, sst et chla - model averaging (fin) -----------------------------------


#--------- relation entre detection chla, sst et effort - model averaging (debut) -----------------------------------

# model selection
fms <- fitList(m1=fm20,m2=fm21,m3=fm22,m4=fm26,m5=fm27,m6=fm28)
modSel(fms)

obscovs %>% head
eff = c(obscovs$effort.km$d.in.grid.1,obscovs$effort.km$d.in.grid.2)
mean_cov <- mean(eff, na.rm = T)
sd_cov <- sd(eff, na.rm = T)
min_cov <- min(eff, na.rm = T)
max_cov <- max(eff, na.rm = T)


mean_cov <- lapply(obscovs,apply, 2, mean, na.rm = T) %>% lapply(.,mean)
sd_cov <- lapply(obscovs,apply, 2, sd, na.rm = T) %>% lapply(.,mean)
min_cov <-lapply(obscovs,min, na.rm = T)
max_cov <- lapply(obscovs,max, na.rm = T)

# model-averaging of detection for a grid of values for effort, with chla and sst to mean
grid_eff_nonst <- seq(min_cov,max_cov,length=1000)
newdata <- data.frame(chla=mean
                      effortkm2=grid_eff_nonst,
)
occ_pred_res <- predict(fms,type='det',newdata=newdata, appendData =TRUE)
occ_pred_eff <- occ_pred_res$Predicted
# effort effect on detection
plot(grid_eff_nonst, occ_pred_eff, xlab= "Effort (km2)", ylab="Probabilit? de d?tection", main="Relation entre d?tection et effort", ylim=c(0,1),las=1,type='n',lwd=3)
polygon(c(rev(grid_eff_nonst), grid_eff_nonst), c(rev(occ_pred_res$upper), occ_pred_res$lower), col = 'grey80', border = NA)
lines(grid_eff_nonst, occ_pred_eff, col="blue",ylim=c(0,1),las=1,type='l',lwd=3)

#--------- relation entre detection et effort - model averaging (fin) -----------------------------------


