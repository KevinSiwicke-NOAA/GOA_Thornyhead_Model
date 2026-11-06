library(afscdata)
library(dplyr)
library(ggplot2)
library(DBI)

# db <- connect() # if you setup keyring this will just work https://afsc-assessments.github.io/afscdata/articles/getting-started.html
db <- DBI::dbConnect (odbc::odbc(),
                      dsn = "akfin",
                      uid = keyring::key_list("akfin")$username,
                      pwd =  keyring::key_get("akfin", keyring::key_list("akfin")$username))
# globals
year <- 2026
species <- 'THDS'
area <- 'goa'

# query data 
abc <- q_specs(year=year, species=species, area=area, db=db, save=F) |> 
  filter(!area_label == "GOA") |> 
  rename(strata = area_label, ABC = acceptable_biological_catch, Year = year) |> 
  mutate(strata = ifelse(strata == "W", "WGOA", 
                         ifelse(strata == "C", "CGOA", "EGOA")))

catch <- dbGetQuery(db, 
                           "select    *
                from      council.comprehensive_blend_ca
                where     species_group_code = 'THDS' and 
                          fmp_area = 'GOA' and year > 1991
                order by  year asc
                ") |> 
  rename_all(tolower) 

DBI::dbDisconnect(db)

area_catch <- catch |> 
  rename(strata = fmp_subarea, Retention = retained_or_discarded) |> 
  mutate(Retention = ifelse(Retention == 'D', 'Discarded', 'Retained'),
         Year = as.numeric(year), 
         strata = ifelse(strata == "WG", "WGOA", 
                         ifelse(strata == "CG", "CGOA", "EGOA"))) |> 
  group_by(Year, strata, Retention) |> 
  summarize(Catch = sum(weight_posted))

catch_disc <- catch |> 
  rename(strata = fmp_subarea, Retention = retained_or_discarded) |> 
  mutate(Retention = ifelse(Retention == 'D', 'Discarded', 'Retained'),
         Year = as.numeric(year), 
         strata = ifelse(strata == "WG", "WGOA", 
                         ifelse(strata == "CG", "CGOA", "EGOA"))) |> 
  group_by(Year, Retention) |> 
  summarize(Catch = sum(weight_posted))

goa_catch <- catch |> 
  mutate(Year = as.numeric(year)) |> 
  group_by(Year) |> 
  summarize(Catch = round(sum(weight_posted), 2))

catch_plot <- left_join(abc, area_catch) |> filter(Year > 2009, Year < 2027)

vals = data.frame()
for(yr in c(2009, 2011, 2013, 2015, 2017, 2018, 2020, 2022, 2024, 2026)) {
  input <- prepare_rema_input(model_name = 'M26',
                              multi_survey = 1,
                              biomass_dat = biomass_dat_3 |> filter(year < yr+1),
                              cpue_dat = cpue_dat |> filter(year < yr+1),
                              sum_cpue_index = TRUE,
                              start_year = 1992,
                              end_year = yr+2,
                              PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1,1,1,1)),
                              q_options = list(
                                pointer_biomass_cpue_strata = c(1, 1,1, 2,2, 2, 3, 3,3),
                                pointer_q_cpue = c(1, 1, 1)),
                              extra_biomass_cv = list(assumption = 'extra_cv'),
                              extra_cpue_cv = list(assumption = 'extra_cv'))
  
  m <- fit_rema(input)
  out <- tidy_rema(m)
  out$parameter_estimates
  
  val <- out$biomass_by_cpue_strata |> filter(year > yr)
  if(yr == 2017) {val <- out$biomass_by_cpue_strata |> filter(year == yr+1)}
  
  vals <- bind_rows(val, vals)
}

ggplot(catch_plot) + 
  geom_col(aes(Year, Catch, fill = Retention), col = 'black') +
  scale_fill_manual(values = c('lightblue', 'steelblue')) +
  geom_point(aes(Year, ABC), col = 'firebrick') +
  geom_line(aes(Year, ABC), col = 'firebrick') +
  scale_y_continuous(expand =c(0,0), limits = c(0, 0.75*0.03*max(vals$pred_uci) + 20)) +
  facet_wrap(~factor(strata, levels = c("WGOA", "CGOA", "EGOA")))  +
  ylab("Catch (t)") +
  theme_bw()

ggsave(filename = paste0(out_path, '/catch_by_area.png'),
       dpi = 300, bg = 'white', units = 'in', height = 6, width = 11)

ggplot(catch_plot) + 
  geom_col(aes(Year, Catch, fill = Retention), col = 'black') +
  scale_fill_manual(values = c('lightblue', 'steelblue')) +
  geom_point(aes(Year, ABC), col = 'firebrick') +
  geom_line(data = vals, aes(year, 0.75*0.03*pred)) +
  geom_point(data = vals, aes(year, 0.75*0.03*pred)) +
  # geom_ribbon(data = vals, aes(year, ymin = 0.75*0.03*pred_lci, ymax = 0.75*0.03*pred_uci), alpha = 0.5) +
  geom_line(aes(Year, ABC), col = 'firebrick') +
  scale_y_continuous(expand =c(0,0), limits = c(0, 0.75*0.03*max(vals$pred_uci) + 20)) +
  facet_wrap(~factor(strata, levels = c("WGOA", "CGOA", "EGOA")))  +
  ylab("Catch (t)") +
  theme_bw()

ggsave(filename = paste0(out_path, '/catch_by_area2.0.png'),
       dpi = 300, bg = 'white', units = 'in', height = 6, width = 11)

ggplot(catch_plot) + 
  geom_col(aes(Year, Catch, fill = Retention)) +
  geom_point(aes(Year, ABC)) +
  geom_line(aes(Year, ABC)) +
  scale_y_continuous(expand =c(0,0), limits = c(0, max(catch_plot$ABC) + 20)) +
  facet_wrap(~factor(strata, levels = c("WGOA", "CGOA", "EGOA")))  +
  ylab("Catch (t)") +
  theme_bw()

ggsave(filename = paste0(out_path, '/catch_by_area.png'),
       dpi = 300, bg = 'white', units = 'in', height = 4, width = 10)

fish_catch <- catch |>
  rename(Gear = agency_gear_code) |> 
  mutate(Year = as.numeric(year)) |> 
  group_by(Year, Gear) |> 
  summarize(Catch = sum(weight_posted))

catch_tab <- catch |>
  rename(Gear = agency_gear_code, Retention = retained_or_discarded) |> 
  mutate(Year = as.numeric(year),
         Gear = ifelse(Gear %in% c('PTR', 'NPT'), 'trawl',
                       ifelse(Gear == 'HAL', 'longline', NA))) |>
  filter(!is.na(Gear)) |> 
  group_by(Year, Gear, Retention) |> 
  summarize(Catch = sum(weight_posted))

ggplot(fish_catch) +
  geom_col(aes(Year, Catch, fill = Gear)) +
  scale_y_continuous(expand =c(0,0), limits = c(0, 1200)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_fishery.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)

target_catch <- catch |>
  rename(Target = trip_target_name) |> 
  mutate(Year = as.numeric(year)) |> 
  group_by(Year, Target) |> 
  summarize(Catch = sum(weight_posted))

target_catch_disc <- catch |>
  filter(retained_or_discarded == 'D') |> 
  rename(Target = trip_target_name) |> 
  mutate(Year = as.numeric(year)) |> 
  group_by(Year, Target) |> 
  summarize(Catch = sum(weight_posted))

ggplot(target_catch) +
  geom_col(aes(Year, Catch, fill = Target)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_target.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)


