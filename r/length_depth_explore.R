library(tidyverse)
library(ggplot2)
library(DBI)
library(keyring)

# I get survey data from AKFIN and my credentials are stored with keyring
db <- "akfin"

channel_akfin <- DBI::dbConnect (odbc::odbc(),
                                 dsn = db,
                                 uid = keyring::key_list(db)$username,
                                 pwd =  keyring::key_get(db, keyring::key_list(db)$username))

lls_len <- dbGetQuery(channel_akfin, 
                      "select    *
                from      afsc.lls_length_summary_view
                where     species_code = '30020' and
                          year > 1991 and country = 'United States'
                order by  year asc") %>% 
  rename_all(tolower) |> 
  filter(stratum > 1) # only 2 fish from 0-100m and some 'unknown' depths at stratum == 0

goa_lls_raw_data <- data.frame(year = rep(lls_len$year, lls_len$frequency),
                               region = rep(lls_len$council_sablefish_mgmt_area, lls_len$frequency),
                               length = rep(lls_len$length, lls_len$frequency),
                               strat = rep(lls_len$stratum_description, lls_len$frequency)) |> 
  filter(!region == 'Bering Sea', !region == 'Aleutians') |> 
  mutate(region = factor(ifelse(region == 'Central Gulf of Alaska', 'CGOA',
                                ifelse(region == 'Western Gulf of Alaska', 'WGOA', 'EGOA')),
                         levels = c('WGOA', 'CGOA', 'EGOA')),
         strat = factor(strat, levels = c("0-100m","101-200m","201-300m","301-400m","401-600m",
                                          "601-800m","801-1000m", "1001-1200m", "1201m +")),
         survey = 'LLS')

ggplot(goa_lls_raw_data, aes(region, length, fill = strat)) +
  geom_boxplot() 

goa_lls_summ <- goa_lls_raw_data |> 
  group_by(year, strat, region) |> 
  summarize(N = n(), mean_len = mean(length))

ggplot(goa_lls_summ, aes(year, mean_len)) +
  geom_point() +
  facet_grid(strat~region)

bts_len <- dbGetQuery(channel_akfin, 
                      "select    *
                from      gap_products.akfin_length_v
                where     species_code = '30020' and
                          SURVEY_DEFINITION_ID = 47 and 
                          year = 1990 and
                          length_mm > 0") %>% 
  rename_all(tolower) |>
  mutate(depth_gear_m = ifelse(is.na(depth_gear_m), depth_m, depth_gear_m),
         depth_gear_m = ifelse(depth_gear_m == 0, depth_m, depth_gear_m))

goa_bts_raw_data <- data.frame(year = rep(bts_len$year, bts_len$frequency),
                               region = rep(bts_len$regulatory_area, bts_len$frequency),
                               length = rep(bts_len$length_mm, bts_len$frequency),
                               depth = rep(bts_len$depth_gear_m, bts_len$frequency)) |> 
  mutate(length = length/10, 
         strat = ifelse(depth < 101, "1-100m",
                        ifelse(depth < 201, "101-200m",
                               ifelse(depth < 301, "201-300m",
                                      ifelse(depth < 401, "301-400m",
                                             ifelse(depth < 601, "401-600m",
                                                    ifelse(depth < 801, "601-800m",
                                                           ifelse(depth < 1001, "801-1000m", ">1000m"))))))),
         region = factor(ifelse(region == 'CENTRAL GOA - NMFS', 'CGOA',
                                ifelse(region == 'CENTRAL GOA - INPFC', 'CGOA',
                                       ifelse(region == 'WESTERN GOA - NMFS', 'WGOA',
                                              ifelse(region == 'WESTERN GOA - INPFC', 'WGOA', 'EGOA')))),
                         levels = c('WGOA', 'CGOA', 'EGOA')),
         survey = 'BTS') |> 
  select(!depth)

ggplot(goa_bts_raw_data, aes(region, length, fill = strat)) +
  geom_boxplot() 

all_lengths <- bind_rows(goa_lls_raw_data, goa_bts_raw_data) |> 
  mutate(strat = factor(strat, levels = c("1-100m","101-200m","201-300m","301-400m","401-600m",
                                   "601-800m","801-1000m", "1001-1200m", "1201m +")))

ggplot(all_lengths, aes(strat, length, fill = survey)) +
  facet_grid(~region) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1)) +
  labs(x = '', y = 'Length (cm)', fill = 'Survey')























lls_cat <- dbGetQuery(channel_akfin, 
                  "select    *
                from      afsc.lls_catch_summary_with_nulls_mv
                where     species_code = '30020' and 
                          year = 2025
                ") %>% 
  rename_all(tolower) |> 
  filter(ineffective < 6) |> 
  mutate(cpue = ifelse(is.na(catch_freq), 0 , catch_freq / (45-ineffective))) |> 
  select(geo = council_sablefish_management_area, cpue, depth = intrpdep) |> 
  mutate(geo = ifelse(geo == 'Central Gulf of Alaska', 'CGOA',
                      ifelse(geo == 'Western Gulf of Alaska', 'WGOA', 'EGOA')))

bins <- seq(0,1100, by =50)
labels <- seq(0,1050, by = 50)
lls_cat$cuts = cut(lls_cat$depth, breaks = bins, labels = labels)

lls_summ <- lls_cat |> 
  group_by(geo, cuts) |> 
  summarize(mean_cpue = mean(cpue))
  
ggplot(lls_summ, aes(cuts, mean_cpue)) + 
  geom_point() +
  facet_wrap(~geo)

# Get bottom trawl survey biomass data
look <- dbGetQuery(channel_akfin, 
                   "select    *
                from      gap_products.akfin_catch_v
                where     species_code = '30020' and 
                          survey_definition_id = 47 and
                          year = 2025
                ") %>% 
  rename_all(tolower) |> 
  select(haul, weight_kg, count, depth = depth_gear_m)

hauls <- dbGetQuery(channel_akfin, 
                   "select    *
                from      gap_products.akfin_haul_v
                where     survey_definition_id = 47 and
                          year = 2025
                ") %>% 
  rename_all(tolower) |> 
  select(haul, depth = depth_gear_m, performance, geo = regulatory_area) 

cat_by_dep <- left_join(hauls, look) |> 
  mutate(weight_kg = ifelse(is.na(weight_kg), 0 , weight_kg),
         count = ifelse(is.na(count), 0, count),
         geo = ifelse(geo == 'CENTRAL GOA - NMFS', 'CGOA',
                      ifelse(geo == 'EASTERN GOA - NMFS', 'EGOA', 'WGOA')))

ggplot(cat_by_dep, aes(depth, count)) + 
  geom_point() +
  facet_wrap(~geo)

ggplot(cat_by_dep, aes(depth, weight_kg)) + 
  geom_point() +
  facet_wrap(~geo)

cat_by_dep$cuts = cut(cat_by_dep$depth, breaks = bins, labels = labels)

bts_summ <- cat_by_dep |> 
  group_by(geo, cuts) |> 
  summarize(mean_cpue = mean(weight_kg))

ggplot(bts_summ, aes(cuts, mean_cpue/max(bts_summ$mean_cpue))) + 
  geom_point() +
  geom_point(data = lls_summ, aes(cuts, mean_cpue/max(lls_summ$mean_cpue)), col = 'red', alpha = 0.5) +
  facet_wrap(~geo)
