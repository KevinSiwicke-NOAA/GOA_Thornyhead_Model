library(dplyr)
library(DBI)
library(keyring)
library(lubridate)
library(ggplot2)
library(akgfmaps)
library(sf)
library(RODBC)

# I get survey data from AKFIN and my credentials are stored with keyring
db <- "akfin"
channel_akfin <- DBI::dbConnect (odbc::odbc(),
                                 dsn = db,
                                 uid = keyring::key_list(db)$username,
                                 pwd =  keyring::key_get(db, keyring::key_list(db)$username))

# Get bottom trawl survey biomass data
bts_cpue <- dbGetQuery(channel_akfin, 
                       "select    *
                from      gap_products.akfin_cpue_v
                where     species_code = 30020 and 
                          survey_definition_id = 47 and
                          year > 2020
                ") %>% 
  rename_all(tolower)

bts_sf <- bts_cpue |>  mutate(lat = (latitude_dd_start + latitude_dd_end)/2, lon = (longitude_dd_start + longitude_dd_end)/2,
                              Region = ifelse(regulatory_area %in% c('WESTERN GOA - INPFC', "WESTERN GOA - NMFS"), 'WGOA',
                                              ifelse(regulatory_area %in% c('EASTERN GOA - INPFC', "EASTERN GOA - NMFS"), 'EGOA', 'CGOA'))) |> 
  group_by(species_code, year, Region, station_number = hauljoin) |> 
  summarize(catch = sum(count), lat = mean(lat), lon = mean(lon)) |> 
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) |> mutate(survey = 'BTS')

# Get longline survey RPWs by station
lls_cat<- dbGetQuery(channel_akfin, 
                     "select    *
                from      afsc.lls_catch_summary_with_nulls_mv
                where     species_code = 30020 and 
                          year in (2021, 2023, 2025)
                ") %>% 
  rename_all(tolower)

lls_sf <- lls_cat |> 
  filter(!council_sablefish_management_area %in% c('Bering Sea', 'Aleutians'), !habitat_type == 'Shallow cross-shelf gully', station_number < 500) |> 
  mutate(catch_freq = ifelse(is.na(catch_freq), 0, catch_freq), lat = (start_latitude + end_latitude)/2, lon = (start_longitude + end_longitude)/2,
         Region = ifelse(council_sablefish_management_area == "Central Gulf of Alaska", 'CGOA',
                         ifelse(council_sablefish_management_area == "Western Gulf of Alaska", 'WGOA', 'EGOA'))) |>
  group_by(species_code, year, Region, station_number) |> 
  summarize(catch = sum(catch_freq), lat = mean(lat), lon = mean(lon)) |> 
  st_as_sf(coords = c('lon', 'lat'), crs = 4326) |> mutate(survey = 'LLS')

goa_trl <- akgfmaps::get_base_layers(select.region = "goa", set.crs = "auto")
goa_ll <- akgfmaps::get_base_layers(select.region = "ll.goa", set.crs = "auto")

cat_sf <- bind_rows(lls_sf, bts_sf)

# gullies in the EGOA are not included in the akgfmaps package, so add them
gully_aa <- read_sf("data/lls_gullies.shp")

lls_geom <- st_union(c(st_geometry(goa_ll$survey.strata), st_geometry(gully_aa)))

# presence/absence
ggplot() +
  geom_sf(data = goa_trl$bathymetry |> filter(DEPTH_M %in% c(200,500,1000)), col = 'grey50', size = 0.2, aes(lty = factor(DEPTH_M))) +
  labs(lty = 'Depth (m)') +
  geom_sf(data = goa_trl$akland) +
  geom_sf(data = goa_trl$survey.strata, col = NA, fill = 'steelblue', alpha = 0.1) +
  geom_sf(data = goa_trl$graticule, color = "grey70", alpha = 0.5) +
  geom_sf(data = lls_geom, col = NA, fill = 'firebrick', alpha = 0.1) +
  geom_sf(data = cat_sf |> filter(catch == 0), aes(fill = Region), pch = 'X', size = 1.5) +
  geom_sf(data = cat_sf |> filter(catch > 0), aes(fill = Region), pch = 21, size = 1, alpha = 0.8) +
  facet_grid(year~survey) +
  coord_sf(xlim = c(-1000000, 1250000), 
           ylim = c(400000, 1200000)) +
  scale_x_continuous(name = "Longitude",
                     breaks = goa_trl$lon.breaks) +
  scale_y_continuous(name = "Latitude",
                     breaks = goa_trl$lat.breaks) +
  scale_fill_viridis_d() +
  theme_bw()

ggsave(paste0(out_path, '/sst_survey_maps.png'), height = 8, width = 12, dpi = 300)

# Observer Catch
zip_url <- "https://apps-afsc.fisheries.noaa.gov/ords/r/fma_ols/files/static/v682/fmap/rockfish.zip"
temp_zip <- tempfile(fileext = ".zip")
temp_dir <- tempdir()
download.file(zip_url, destfile = temp_zip, mode = "wb")
unzipped_files <- unzip(temp_zip, exdir = temp_dir)
data <- read.csv(unzipped_files) |> 
  rename(lat = LAT400SQKM, lon = LON400SQKM) |> 
  rename_all(tolower) |> 
  filter(species == 'Shortspine Thornyhead', year > 2022) |> 
  mutate(gear = ifelse(gear %in% c('POT', 'LGL'), 'Fixed', 'Trawl')) |> 
  group_by(year, lat, lon, gear) |> 
  summarize(weight = sum(kg)/ 1000) |> 
  ungroup()

fsh_cat_sf <- data |> st_as_sf(coords = c('lon', 'lat'), crs = 4326) |> 
  st_transform(st_crs(goa_trl$survey.strata)) 

goa_points <- st_filter(fsh_cat_sf, goa_trl$survey.strata, .predicate = st_intersects)
cellsize = 20000
half_size = 10000

cells <- st_as_sfc(lapply(st_geometry(goa_points), function(pt) {
  coords <- st_coordinates(pt)
  xmin <- coords[1] - half_size
  ymin <- coords[2] - half_size
  xmax <- coords[1] + half_size
  ymax <- coords[2] + half_size
  st_polygon(list(matrix(c(
    xmin, ymin,
    xmax, ymin,
    xmax, ymax,
    xmin, ymax,
    xmin, ymin
  ), ncol = 2, byrow = TRUE)))
}))

sf_cells <- st_sf(goa_points[setdiff(names(goa_points), "geometry")], geometry = cells, crs = st_crs(goa_points))

ggplot() +
  geom_sf(data = goa_trl$bathymetry |> filter(DEPTH_M %in% c(200,500,1000)), col = 'grey50', size = 0.2, aes(lty = factor(DEPTH_M))) +
  labs(lty = 'Depth (m)', fill = 'Catch (t)') +
  geom_sf(data = goa_trl$akland) +
  geom_sf(data = goa_trl$survey.strata, col = NA, fill = 'steelblue', alpha = 0.1) +
  geom_sf(data = goa_trl$graticule, color = "grey70", alpha = 0.5) +
  geom_sf(data = lls_geom, col = NA, fill = 'firebrick', alpha = 0.1) +
  geom_sf(data = sf_cells, aes(fill = weight), pch = 21, alpha = 0.5) +
  facet_grid(year~gear) +
  coord_sf(xlim = c(-1100000, 1300000), 
           ylim = c(320000, 1200000)) +
  scale_x_continuous(name = "Longitude",
                     breaks = goa_trl$lon.breaks) +
  scale_y_continuous(name = "Latitude",
                     breaks = goa_trl$lat.breaks) +
  scale_fill_viridis_c() +
  theme_bw()

ggsave(paste0(out_path, '/sst_observer_maps.png'), height = 8, width = 12, dpi = 300)
