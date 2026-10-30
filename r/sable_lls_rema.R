# install.packages("devtools")
# devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)

library(tidyverse)
library(DBI)
library(keyring)
library(DBI)
library(corrplot)
library(rema)

# I get survey data from AKFIN and my credentials are stored with keyring
db <- "akfin"

channel_akfin <- dbConnect (odbc::odbc(),
                                 dsn = db,
                                 uid = keyring::key_list(db)$username,
                                 pwd =  keyring::key_get(db, keyring::key_list(db)$username))

# Get longline survey RPWs, which is an index for biomass. This will approximate distribution of biomass across management areas.
rpw <- dbGetQuery(channel_akfin, 
                  "select    *
                from      afsc.lls_area_rpn_3_to_7_depred
                where     country = 'United States' and 
                          exploitable = 1
                order by  year asc
                ") |> 
  rename_all(tolower) |> 
  filter(!geographic_area_name %in% c('NW Aleutians slope', 'SW Aleutians slope')) # Remove the WAI

# Calling the index 'biomass' keeps rema from yelling at you, as it expects biomass from the trawl survey or something like that.
lls_index <- rpw |> 
  filter(year > 1989, !(year == 1996 & council_sablefish_management_area == 'Aleutians')) |> # Remove AI 96, which is incomplete in current database. There are fixed values to get this back correctly, but not done at this time.
  group_by(year, strata = council_sablefish_management_area) |> 
  mutate(strata = ifelse(strata == 'Western Gulf of Alaska', 'WGOA',
                         ifelse(strata == 'Central Gulf of Alaska', 'CGOA',
                                ifelse(strata == 'West Yakutat', 'WY', 
                                       ifelse(strata == 'East Yakutat/Southeast', 'EYSE', 
                                              ifelse(strata == 'Aleutians', 'AI', 
                                                     ifelse(strata == 'Bering Sea', 'BS', NA))))))) |> 
  summarize(biomass = sum(rpw, na.rm = TRUE),
            cv = sqrt(sum(rpw_var, na.rm = TRUE)) / biomass)

DBI::dbDisconnect(channel_akfin)

# build data for model
lls_dat <- left_join(data.frame('year' = rep(unique(lls_index$year), each = 6), 'strata' = rep(c('AI', 'BS', 'WGOA', 'CGOA', 'WY', 'EYSE'), length(unique(lls_index$year)))), lls_index, by = c('year', 'strata')) 

input <- prepare_rema_input(model_name = 'lls_rema',
                            biomass_dat = lls_dat, # to edit the input years of data, just add |> filter(year < YYYY)
                            start_year = 1990,
                            end_year = 2027, # change as needed depending how far past the end of the filtered data you want to estimate
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1))) # All regions share a single process error

m_lls <- fit_rema(input)
out <- tidy_rema(m_lls)
out$parameter_estimates

# A few canned plots from package, but can customize if desired
rema_plots <- plot_rema(out) 

rema_plots$biomass_by_strata + theme_bw()
rema_plots$total_predicted_biomass + theme_bw()
rema_plots$proportion_biomass_by_strata + theme_bw()

# apportionment ----
appo <- out$biomass_by_strata |>
  group_by(strata, year) |>
  summarize(stratum_biomass = sum(pred)) |>
  group_by(year) |>
  mutate(total_biomass = sum(stratum_biomass)) |>
  ungroup() |>
  mutate(proportion = stratum_biomass / total_biomass,
         region = ifelse(strata %in% c('CGOA', 'WGOA', 'WY', 'EYSE'), 'GOA', strata))

# Plotting out results for apportionment
# All together
ggplot(appo, aes(year, stratum_biomass, fill = strata)) + 
  geom_col() + 
  labs(x = 'Year', y = 'RPWs', fill = 'Area') +
  scale_y_continuous(limits = c(0, max(appo$total_biomass)+100000), expand = c(0,0)) +
  theme_bw()

# By FMP
ggplot(appo, aes(year, stratum_biomass, fill = strata)) + 
  geom_col() + 
  labs(x = 'Year', y = 'RPWs', fill = 'Area') +
  facet_wrap(~region) +
  scale_y_continuous(limits = c(0, 1800000), expand = c(0,0)) +
  theme_bw()

# By management area
ggplot(appo, aes(year, stratum_biomass)) + 
  geom_col() + 
  labs(x = 'Year', y = 'RPWs') +
  facet_wrap(~strata) + 
  scale_y_continuous(limits = c(0, max(appo$stratum_biomass)+50000), expand = c(0,0)) +
  theme_bw()

# Same plots, but instead of in RPWs, as proportions
ggplot(appo, aes(year, proportion, fill = strata)) + 
  geom_col() + 
  labs(x = 'Year', y = 'Proportion', fill = 'Area') +
  theme_bw()

ggplot(appo, aes(year, proportion, fill = strata)) + 
  geom_col() +
  labs(x = 'Year', y = 'Proportion', fill = 'Area') +
  facet_wrap(~region) + 
  theme_bw()

ggplot(appo, aes(year, proportion)) + 
  geom_col() + 
  labs(x = 'Year', y = 'Proportion') +
  facet_wrap(~strata) + 
  theme_bw()

# Just to get a sense of correlation by areas
for_cor <- out$biomass_by_strata |> 
  select(strata, year, pred) |>
  pivot_wider(names_from = strata, values_from = pred) |> 
  select(-year)

cor_matrix <- cor(for_cor)

corrplot(cor_matrix, method = "color", type = "upper", addCoef.col = "black")
