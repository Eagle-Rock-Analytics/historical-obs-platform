## List of buoys with wind data

library(rnoaa)
library(rvest)
library(tidyverse)
library(tmap)
library(sf)

## All buoys
x <- buoy_stations()

## Off coast of CA
## east of -137 W
## within 44N-29N
ca_buoys <- x %>%
  filter(lat > 29 & lat < 45) %>%
  filter(lon > -137 & lon < -114)

## Get continuous winds and peak winds buoys
cwinds <- buoys('cwind')

pwinds <- buoys('pwind')

## Stations in continuous winds and peak winds dataset in geographic filter
cwinds_ca <- cwinds %>%
  filter(id %in% ca_buoys$station)

pwinds_ca <- pwinds %>%
  filter(id %in% ca_buoys$station)
## These are the same 25

## Find what years data is available for each buoy 
cwinds_years <- cwinds_ca %>%
  mutate(years_html = map(url, ~{
    web <- read_html(.)
    
    yrs <- web %>% 
      html_table(fill = T) %>%
      data.frame() %>%
      mutate(year1 = word(Dataset, 2, 2, sep = "c"),
             year = word(year1, 1, 1, sep = "\\.")) %>%
      filter(year != "ml", !is.na(year))
    
    yrs$year
  }))

## Unnest

cwinds_unnest <- cwinds_years %>%
  unnest(c("years_html"))

stations_per_year <- cwinds_unnest %>%
  group_by(years_html) %>%
  summarize(n_stations = n_distinct(id)) %>%
  filter(years_html != "") %>%
  rename(year = "years_html")
write.csv(stations_per_year, "buoy_wind_data/buoys_w_wind_per_year.csv", row.names =F)

## Map of US states
states <- read_sf("map_data/cb_2018_us_state_5m.shp")
west_coast <- states %>%
  filter(NAME %in% c("California", "Oregon", "Washington", "Nevada", "Idaho"))

## buoy points, number of years of data
buoys_sf <- ca_buoys %>%
  right_join(cwinds_unnest, by = c("station" = "id")) %>%
  filter(years_html != "") %>%
  group_by(station, lat, lon) %>%
  summarize(n_years = n_distinct(years_html)) %>%
  st_as_sf(coords = c("lon", "lat")) %>%
  st_set_crs(st_crs(states))

## Map of buoy sites

crop_west <- states %>%
  st_crop(c(xmin = -140, ymin = 29, ymax = 45, xmax = -110))

buoy_map <- tm_shape(buoys_sf) + tm_symbols(col = "n_years", size = 0.5, palette = "Blues", border.col = "black", title.col = "Years of data") +
  tm_shape(crop_west) + tm_polygons()
tmap_save(buoy_map, "buoy_wind_data/map_of_buoy.pdf")
