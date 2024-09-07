library(tidyverse)
library(DBI)
library(keyring)

# I get survey data from AKFIN and my credentials are stored with keyring
db <- "akfin"

channel_akfin <- DBI::dbConnect (odbc::odbc(),
                                 dsn = db,
                                 uid = keyring::key_list(db)$username,
                                 pwd =  keyring::key_get(db, keyring::key_list(db)$username))


# Get longline survey RPWs
rpw <- dbGetQuery(channel_akfin, 
                "select    *
                from      afsc.lls_area_rpn_all_strata
                where     species_code = '30020' and 
                          fmp_management_area = 'GOA' and 
                          exploitable = 1 and 
                          country = 'United States'
                order by  year asc
                ") %>% 
  rename_all(tolower)

cpue <- rpw %>% 
  filter(year > 1991) %>% 
  group_by(year, strata = council_management_area) %>% 
  mutate(strata = ifelse(strata == 'Western Gulf of Alaska', 'WGOA',
                         ifelse(strata == 'Central Gulf of Alaska', 'CGOA',
                                ifelse(strata == 'Eastern Gulf of Alaska', 'EGOA', NA)))) %>% 
  summarize(cpue = sum(rpw, na.rm = TRUE),
            cv = sqrt(sum(rpw_var, na.rm = TRUE)) / cpue)

cpue_dat <- left_join(data.frame('year' = rep(unique(cpue$year), each = 3), 'strata' = rep(c('WGOA', 'CGOA', 'EGOA'), length(unique(cpue$year)))), cpue, by = c('year', 'strata')) 

# ggplot(cpue_dat, aes(year, cpue)) +
#    geom_line() +
#    facet_wrap(~strata)

# Get bottom trawl survey biomass data
biom <- dbGetQuery(channel_akfin, 
                       "select    *
                from      gap_products.akfin_biomass
                where     species_code = '30020' and 
                          survey_definition_id = 47
                order by  year asc
                ") %>% 
  rename_all(tolower) |> 
  mutate(biomass_var = ifelse(is.na(biomass_var), (0.5 * biomass_mt) ^ 2, 
                              ifelse(biomass_var == 0 & biomass_mt > 0, (0.5 * biomass_mt) ^ 2, biomass_var)))

biomass <- biom %>% 
  mutate(strata = ifelse(area_id %in% c(10:13, 110:112, 210, 310), 'WGOA (0-500 m)',
                         ifelse(area_id %in% c(20:35, 120:134, 220:232, 32, 320, 330), 'CGOA (0-500 m)',
                                ifelse(area_id %in% c(40:50, 140:151, 240:251, 340:351), 'EGOA (0-500 m)',
                                       ifelse(area_id == 410, 'WGOA (501-700 m)',
                                              ifelse(area_id == 510, 'WGOA (701-1000 m)',
                                                     ifelse(area_id %in% c(420, 430), 'CGOA (501-700 m)',
                                                            ifelse(area_id %in% c(520, 530), 'CGOA (701-1000 m)',
                                                                   ifelse(area_id %in% c(440, 450), 'EGOA (501-700 m)',
                                                                          ifelse(area_id %in% c(540, 550), 'EGOA (701-1000 m)', NA)))))))))) %>% 
  filter(!is.na(strata)) %>% 
  group_by(year, strata) %>% 
  summarize(n = sum(n_haul), biomass = sum(biomass_mt, na.rm = TRUE),
            cv = sqrt(sum(biomass_var, na.rm = TRUE))/biomass) 

biomass_dat <- left_join(data.frame('year' = rep(unique(biomass$year), each = 9), 'strata' = rep(unique(biomass$strata), length(unique(biomass$year)))), biomass, by = c('year', 'strata')) %>% 
  mutate(cv = ifelse (cv == 0, 0.5, cv))  

model_yrs <- 1990:YEAR

# ggplot(biomass_dat, aes(year, biomass)) +
#   geom_line() +
#   facet_wrap(~strata)

# This is the data that is brought into rema
model_dat <- list('biomass_dat' = biomass_dat, 'cpue_dat' = cpue_dat, 
                  'model_yrs' = model_yrs)

DBI::dbDisconnect(channel_akfin)
