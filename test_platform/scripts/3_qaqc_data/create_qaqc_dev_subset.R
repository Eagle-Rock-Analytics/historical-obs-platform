## Eval training dataset for QA/QC

library(tidyverse)
library(sf)
library(raster)
library(lhs)

# Read in list of cleaned stations
stns <- read_csv('cleaned_stations.csv')

# Read in bioclim vars
bioclim_files <- list.files('map_data/')[grepl('bio', list.files('map_data/'))]

bioclim <- stack(paste0('map_data/',bioclim_files))

# Extract bioclim at stations
stns_sf <- stns %>%
  st_as_sf(coords = c('longitude', 'latitude')) %>%
  st_set_crs(4236)

bioclim_stns <- extract(bioclim, stns_sf)

# Join bioclim vars to stns df
stns_allvars <- bind_cols(stns_sf, bioclim_stns)

# Keep bio5,6,13,14 (max temp of warmest month, min temp of coldest month, 
# precip driest month, precip wettest month)
stns_subset <- stns_allvars %>%
  dplyr::select(name, elevation, network, 
                wc2.1_30s_bio_5, wc2.1_30s_bio_6, wc2.1_30s_bio_13, wc2.1_30s_bio_14)

# 250 stations across 1 dimension, maximin criteria LHS
lhs_samp <- maximinLHS(250, 1)

# function: convert lhs index to a new variable domain
# inputs: LHS output value to convert
# inputs: vector with values of the new variable
transform_lhs <- function(lhs_index, variable_vector) {
  min_var <- min(variable_vector, na.rm = T)
  max_var <- max(variable_vector, na.rm = T)
  range_var <- max_var - min_var
  
  conv <- min_var + lhs_index*range_var
  
  return(conv)
}

# match lhs positions to stations for each variable
elev <- pmap(list(lhs_samp), 
     ~which.min(abs(stns_subset$elevation - transform_lhs(., stns_subset$elevation)))) %>%
  unlist()

bio5 <- pmap(list(lhs_samp), 
             ~which.min(abs(stns_subset$wc2.1_30s_bio_5 - transform_lhs(., stns_subset$wc2.1_30s_bio_5)))) %>%
  unlist()

bio6 <- pmap(list(lhs_samp), 
             ~which.min(abs(stns_subset$wc2.1_30s_bio_6 - transform_lhs(., stns_subset$wc2.1_30s_bio_6)))) %>%
  unlist()

bio13 <- pmap(list(lhs_samp), 
              ~which.min(abs(stns_subset$wc2.1_30s_bio_13 - transform_lhs(., stns_subset$wc2.1_30s_bio_13)))) %>%
  unlist()

bio14 <- pmap(list(lhs_samp), 
              ~which.min(abs(stns_subset$wc2.1_30s_bio_14 - transform_lhs(., stns_subset$wc2.1_30s_bio_14)))) %>%
  unlist()

# unique list of stations
stns_lhs <- unique(c(elev, bio5, bio6, bio13, bio14))

# get stations in final list
stns_final <- stns_subset[stns_lhs, ]

## Create plots showing full distribution of key variables and testing subset distributions
theme_set(theme_classic())

elev <- ggplot(stns_subset, aes(x = elevation)) + 
  geom_density() + ylab('density') +
  geom_point(data = stns_final, aes(x = elevation, y = 0), col = 'skyblue', alpha = 0.5)

bio5 <- ggplot(stns_subset, aes(x = wc2.1_30s_bio_5)) + 
  geom_density() + ylab('density') + xlab('maximum temp warmest month') +
  geom_point(data = stns_final, aes(x = wc2.1_30s_bio_5, y = 0), col = 'skyblue', alpha = 0.5)

bio6 <- ggplot(stns_subset, aes(x = wc2.1_30s_bio_6)) + 
  geom_density() + ylab('density') + xlab('minimum temp coldest month') +
  geom_point(data = stns_final, aes(x = wc2.1_30s_bio_6, y = 0), col = 'skyblue', alpha = 0.5)

bio13 <- ggplot(stns_subset, aes(x = wc2.1_30s_bio_13)) + 
  geom_density() + ylab('density') + xlab('precip driest month') +
  geom_point(data = stns_final, aes(x = wc2.1_30s_bio_13, y = 0), col = 'skyblue', alpha = 0.5)

bio14 <- ggplot(stns_subset, aes(x = wc2.1_30s_bio_14)) + 
  geom_density() + ylab('density') + xlab('precip wettest month') +
  geom_point(data = stns_final, aes(x = wc2.1_30s_bio_14, y = 0), col = 'skyblue', alpha = 0.5)

cowplot::plot_grid(elev, bio5, bio6, bio13, bio14, nrow = 3)

## Write df with final station set
coords <- st_coordinates(stns_final) %>%
  data.frame() %>%
  rename('longitude' = X, 'latitude' = Y)

stns_to_write <- stns_final %>%
  st_set_geometry(NULL) %>%
  dplyr::select(name, elevation, network) %>%
  bind_cols(coords)

write.csv(stns_to_write, "test_platform/scripts/3_qaqc_data/qaqc_training_station_list.csv", row.names = F)
