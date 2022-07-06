library(httr)
library(sf)

# Pull marine shapefile
url <- parse_url("https://services2.arcgis.com/z5VWEXupJpBJnZRq/ArcGIS/rest/services")
url$path <- paste(url$path, "WECC_Informational_MarineCoastal_Boundary/FeatureServer/0/query", sep = "/")
url$query <- list(where = "in_WECC = 'Y'",
                  outFields = "*",
                  returnGeometry = "true",
                  f = "geojson")
request <- build_url(url)

marine_wecc <- st_read(request)

write_sf(marine_wecc, "test_platform/data/maps/WECC_Informational_MarineCoastal_Boundary_marine.shp")

# Pull terrestrial shapefile
url <- parse_url("https://services2.arcgis.com/z5VWEXupJpBJnZRq/ArcGIS/rest/services")
url$path <- paste(url$path, "WECC_Informational_MarineCoastal_Boundary/FeatureServer/1/query", sep = "/")
url$query <- list(where = "in_WECC = 'Y'",
                  outFields = "*",
                  returnGeometry = "true",
                  f = "geojson")
request <- build_url(url)

land_wecc <- st_read(request)

write_sf(land_wecc, "test_platform/data/maps/WECC_Informational_MarineCoastal_Boundary_land.shp")
