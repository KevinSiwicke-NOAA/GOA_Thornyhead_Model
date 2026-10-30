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
rpns <- dbGetQuery(channel_akfin, 
                  "select    *
                from      afsc.lls_area_rpn_all_strata
                where     species_code in (30030, 30020) and 
                          fmp_management_area = 'GOA' and 
                          exploitable = 1 and 
                          country = 'United States'
                order by  year asc
                ") |> 
  rename_all(tolower) |> 
  mutate(sp = ifelse(species_code == 30020, 'SST', 'LST')) |> 
  group_by(year, sp) |> 
  summarize(rpn = sum(rpn, na.rm = T)) |> 
  filter(year > 1991) |> 
  pivot_wider(names_from = sp, values_from = rpn)

ggplot(rpns, aes(year, LST/SST)) + geom_point() + theme_bw()

biom <- dbGetQuery(channel_akfin, 
                   "select    *
                from      gap_products.akfin_biomass_v
                where     species_code in (30030, 30020) and 
                          area_id = 99903
                order by  year asc
                ") |> 
  rename_all(tolower) |>
  mutate(sp = ifelse(species_code == 30020, 'SST', 'LST')) |> 
  filter(year > 1989) |> 
  select(year, sp, biomass_mt) |> 
  pivot_wider(names_from = sp, values_from = biomass_mt)

DBI::dbDisconnect(channel_akfin)

ggplot(biom, aes(year, LST/SST)) + geom_point() + theme_bw()
