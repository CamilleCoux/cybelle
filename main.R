# script based on dauphin_occ2

library(magrittr)
library(tidyr)
library(dplyr)
library(AICcmodavg)
library(unmarked)
library(ggplot2)
library(sf)

theme_set(theme_minimal())

# import obs data
track <- read.csv2("/Users/camillecoux/Dropbox/cybelle/data/track_2009_2020_simple.csv", row.names = 1)

# import grid:
m_grid <-  st_read("data/env_data/mediterranee shapefiles/m_grid.shp", crs=4326)




# select observations for the striped dolphin species, i.e. "dauphn bleu et blanc" in french:
presences <- track[grep("Dauphin Bleu", track$species),]
# selection of all absences  
absences <- track[which(is.na(track$species)),]
# write.csv2(absences, "../cartes/absences.csv ")

# merge
striped_dolphin <- rbind(presences, absences)
# vire les obs ponctuelles
striped_dolphin <- striped_dolphin[-grep("ponctuelle", striped_dolphin$protocole),]

# keep observations from the summer months only
presences$month %>% table %>% barplot # donc on prend les 4 mois d'été
striped_dolphin <- striped_dolphin %>%
  filter(month %in% 6:9)

# one model per year : keep data from 2015 to 2020
table(striped_dolphin$month, striped_dolphin$year)
striped_dolphin <- striped_dolphin %>% filter(year %in% c(2015:2020))


