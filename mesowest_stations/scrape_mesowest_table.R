### Scrape table of station networks from Meso West
#https://mesowest.utah.edu/cgi-bin/droman/raws_ca_monitor.cgi?state=CA&rawsflag=290&timeobs=12&orderby=n&type=0&refrsh=1&stnorder=0

library(tidyverse)
library(rvest)
library(janitor)

content <- read_html("mesowest_stations/MesoWest Weather Summary.html")

table_info <- content %>% html_table(fill = T)

wx_table <- table_info[[6]]

wx_cleannames <- wx_table[-c(1:3), ] %>%
  clean_names()

colnames(wx_cleannames) <- wx_cleannames[1,]

test <- wx_cleannames %>%
  filter(grepl("County", Elev))

test <- wx_cleannames %>%
  filter(grepl("^$", Elev))

wx_cleanrows <- wx_cleannames[-1,] %>%
  filter(!grepl("County", Elev), # RM rows that are county headers for subtables on website
         !grepl("Elev", Elev), # RM rows that are col header repeats for subtables
         !grepl("Time", LOCAL), # RM unit rows repeats for subtables
         !grepl("^$", LOCAL)) # RM rows that are all empty strings separating subtables

networks <- wx_cleanrows %>%
  group_by(Network) %>% 
  summarize(stations = n_distinct(Station)) %>%
  arrange(desc(stations))

write.csv(networks, "mesowest_stations/networks_n_stations.csv", row.names = F)
