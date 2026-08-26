library(afscdata)
library(dplyr)
library(ggplot2)

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
                          fmp_area = 'GOA' and year > 2009
                order by  year asc
                ") |> 
  rename_all(tolower) 

DBI::dbDisconnect(db)

area_catch <- catch |> 
  rename(strata = fmp_subarea, Retention = retained_or_discarded) |> 
  mutate(Year = as.numeric(year), strata = ifelse(strata == "WG", "WGOA", 
                                                  ifelse(strata == "CG", "CGOA", "EGOA"))) |> 
  group_by(Year, strata, Retention) |> 
  summarize(Catch = sum(weight_posted))

catch_plot <- left_join(abc, area_catch) |> filter(Year > 2009, Year < 2027)

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

ggplot(fish_catch) +
  geom_col(aes(Year, Catch, fill = Gear)) +
  scale_y_continuous(expand =c(0,0), limits = c(0, 1200)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_fishery.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)

target_catch <- catch |>
  rename(Target = trip_target_code) |> 
  mutate(Year = as.numeric(year)) |> 
  group_by(Year, Target) |> 
  summarize(Catch = sum(weight_posted))

ggplot(target_catch) +
  geom_col(aes(Year, Catch, fill = Target)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_target.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)


