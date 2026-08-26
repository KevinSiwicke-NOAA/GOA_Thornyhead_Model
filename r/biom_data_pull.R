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
                ") |> 
  rename_all(tolower)

cpue <- rpw |> 
  filter(year > 1991) |> 
  group_by(year, strata = council_management_area) |> 
  mutate(strata = ifelse(strata == 'Western Gulf of Alaska', 'WGOA',
                         ifelse(strata == 'Central Gulf of Alaska', 'CGOA',
                                ifelse(strata == 'Eastern Gulf of Alaska', 'EGOA', NA)))) |> 
  summarize(cpue = sum(rpw, na.rm = TRUE),
            cv = sqrt(sum(rpw_var, na.rm = TRUE)) / cpue)

cpue_dat <- left_join(data.frame('year' = rep(unique(cpue$year), each = 3), 'strata' = rep(c('WGOA', 'CGOA', 'EGOA'), length(unique(cpue$year)))), cpue, by = c('year', 'strata')) 

# ggplot(cpue_dat, aes(year, cpue)) +
#    geom_line() +
#    facet_wrap(~strata)

# Get bottom trawl survey biomass data
biom <- dbGetQuery(channel_akfin, 
                       "select    *
                from      gap_products.akfin_biomass_v
                where     species_code = '30020' and 
                          survey_definition_id = 47
                order by  year asc
                ") |> 
  rename_all(tolower)

biomass_3 <- biom |> 
  mutate(biomass_var = ifelse(is.na(biomass_var), (0.5 * biomass_mt) ^ 2, 
                              ifelse(biomass_var == 0 & biomass_mt > 0, (0.5 * biomass_mt) ^ 2, biomass_var)),
         strata = ifelse(area_id %in% c(10:13, 110:112, 210, 310), 'WGOA_0_500',
                         ifelse(area_id %in% c(20:35, 120:134, 220:232, 32, 320, 330), 'CGOA_0_500',
                                ifelse(area_id %in% c(40:50, 140:151, 240:251, 340:351), 'EGOA_0_500',
                                       ifelse(area_id == 410, 'WGOA_501_700',
                                              ifelse(area_id == 510, 'WGOA_701_1000',
                                                     ifelse(area_id %in% c(420, 430), 'CGOA_501_700',
                                                            ifelse(area_id %in% c(520, 530), 'CGOA_701_1000',
                                                                   ifelse(area_id %in% c(440, 450), 'EGOA_501_700',
                                                                          ifelse(area_id %in% c(540, 550), 'EGOA_701_1000', NA)))))))))) |> 
  filter(!is.na(strata)) |> 
  group_by(year, strata) |> 
  summarize(n = sum(n_haul), biomass = sum(biomass_mt, na.rm = TRUE),
            cv = sqrt(sum(biomass_var, na.rm = TRUE))/biomass) 

# Includes 2025 redesign that can no longer be separated along 0-500 and 501-700, so combine to 0-700.
biomass_2 <- biom |> 
  mutate(biomass_var = ifelse(is.na(biomass_var), (0.5 * biomass_mt) ^ 2, 
                              ifelse(biomass_var == 0 & biomass_mt > 0, (0.5 * biomass_mt) ^ 2, biomass_var)),
         strata = ifelse(area_id %in% c(10:15, 110:113, 210:211, 310, 410), 'WGOA_0_700',
                         ifelse(area_id %in% c(20:38, 120:136, 220:232, 32, 320:321, 330, 420, 430), 'CGOA_0_700',
                                ifelse(area_id %in% c(40:51, 140:152, 240:252, 340:352, 440, 450), 'EGOA_0_700',
                                       ifelse(area_id %in% c(510:511), 'WGOA_701_1000',
                                              ifelse(area_id %in% c(520:521, 530:531), 'CGOA_701_1000',
                                                     ifelse(area_id %in% c(540:541, 550:551), 'EGOA_701_1000', NA))))))) |> 
  filter(!is.na(strata)) |> 
  group_by(year, strata) |> 
  summarize(n = sum(n_haul), biomass = sum(biomass_mt, na.rm = TRUE),
            cv = sqrt(sum(biomass_var, na.rm = TRUE))/biomass) 

biomass_dat_3strat <- left_join(data.frame('year' = rep(unique(biomass_3$year), each = 9), 'strata' = rep(unique(biomass_3$strata), length(unique(biomass_3$year)))), biomass_3, by = c('year', 'strata')) |> 
  mutate(cv = ifelse (cv == 0, 0.5, cv))  

mns_5_7 <- biomass_dat_3strat |> group_by(strata) |> summarize(strat_mn = mean(biomass, na.rm = T))
w_5_7 <- mns_5_7 |> filter(strata == 'WGOA_501_700') |> pull(strat_mn)
c_5_7 <- mns_5_7 |> filter(strata == 'CGOA_501_700') |> pull(strat_mn)
e_5_7 <- mns_5_7 |> filter(strata == 'EGOA_501_700') |> pull(strat_mn)

biomass_3_wide <- biomass_dat_3strat |> 
  select(-c(n, cv)) |> 
  pivot_wider(names_from = strata, values_from = biomass) |> 
  mutate(WGOA_501_700 = ifelse(is.na(WGOA_501_700), w_5_7, WGOA_501_700),
         CGOA_501_700 = ifelse(is.na(CGOA_501_700), c_5_7, CGOA_501_700),
         EGOA_501_700 = ifelse(is.na(EGOA_501_700) & !is.na(EGOA_0_500), e_5_7, EGOA_501_700))
# 
# biomass_3_wide <- biomass_dat_3strat |> 
#   select(-c(n, cv)) |> 
#   pivot_wider(names_from = strata, values_from = biomass) |> 
#   mutate(wgoa_per_5_7 = WGOA_501_700 / WGOA_0_500,
#          cgoa_per_5_7 = CGOA_501_700 / CGOA_0_500,
#          egoa_per_5_7 = EGOA_501_700 / EGOA_0_500,
#          WGOA_501_700 = ifelse(is.na(WGOA_501_700), mean(wgoa_per_5_7, na.rm = T) * WGOA_0_500, WGOA_501_700),
#          CGOA_501_700 = ifelse(is.na(CGOA_501_700), mean(cgoa_per_5_7, na.rm = T) * CGOA_0_500, CGOA_501_700),
#          EGOA_501_700 = ifelse(is.na(EGOA_501_700), mean(egoa_per_5_7, na.rm = T) * EGOA_0_500, EGOA_501_700)) |> 
#   select(!c(wgoa_per_5_7, cgoa_per_5_7, egoa_per_5_7))

biomass_3strat_impute <- biomass_3_wide |> 
  pivot_longer(cols = contains("GOA"), names_to = 'strata', values_to = 'biomass') |> 
  left_join(biomass_dat_3strat |> select(year, strata, cv)) |> 
  mutate(cv = ifelse(!is.na(biomass) & is.na(cv), 0.5, cv)) |> 
  filter(year < 2025) # 2025 doesn't cooperate because of survey changes, so just add it from the non-impute version

combine_strata_pair <- function(pair, data) {
  # Filter data for just the two strata
  d_pair <- data %>% filter(strata %in% pair)
  
  # Calculate combined metrics
  d_pair %>%
    mutate(var = (biomass * cv)^2) |> 
    group_by(year) |> 
    summarise(
      biomass = sum(biomass),
      cv = sqrt(sum(var)) / biomass
    ) |> 
    ungroup()
}

w_2strat_impute <- combine_strata_pair(pair = c("WGOA_0_500", "WGOA_501_700"), biomass_3strat_impute) |> mutate(strata = 'WGOA_0_700') |> select(year, strata, biomass, cv)
c_2strat_impute <- combine_strata_pair(pair = c("CGOA_0_500", "CGOA_501_700"), biomass_3strat_impute) |> mutate(strata = 'CGOA_0_700') |> select(year, strata, biomass, cv)
e_2strat_impute <- combine_strata_pair(pair = c("EGOA_0_500", "EGOA_501_700"), biomass_3strat_impute) |> mutate(strata = 'EGOA_0_700') |> select(year, strata, biomass, cv)
  
biomass_2strat_impute <- bind_rows(w_2strat_impute, c_2strat_impute, e_2strat_impute, biomass_3strat_impute |> filter(stringr::str_detect(strata, "701_1000"))) |> 
  filter(year < 2025) # 2025 doesn't cooperate because of survey changes, so just add it from the non-impute version

biomass_dat_2strat <- left_join(data.frame('year' = rep(unique(biomass_2$year), each = 6), 'strata' = rep(unique(biomass_2$strata), length(unique(biomass_2$year)))), biomass_2, by = c('year', 'strata')) |> 
  mutate(cv = ifelse (cv == 0, 0.5, cv))  

biomass_2strat_impute <- bind_rows(biomass_2strat_impute, biomass_dat_2strat |> filter(year == 2025) )

model_yrs <- 1990:YEAR

# ggplot(biomass_dat, aes(year, biomass)) +
#   geom_line() +
#   facet_wrap(~strata)

# This is the data that is brought into rema
model_dat <- list('biomass_dat_3strat' = biomass_dat_3strat, 
                  'biomass_dat_3strat_impute' = biomass_3strat_impute,
                  'biomass_dat_2strat' = biomass_dat_2strat, 
                  'biomass_dat_2strat_impute' = biomass_2strat_impute, 
                  'cpue_dat' = cpue_dat, 
                  'model_yrs' = model_yrs)

DBI::dbDisconnect(channel_akfin)
