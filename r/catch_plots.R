library(afscdata)
library(dplyr)
library(ggplot2)

# db <- connect() # if you setup keyring this will just work https://afsc-assessments.github.io/afscdata/articles/getting-started.html
db <- DBI::dbConnect (odbc::odbc(),
                      dsn = "akfin",
                      uid = keyring::key_list("akfin")$username,
                      pwd =  keyring::key_get("akfin", keyring::key_list("akfin")$username))
# globals
year <- 2024
species <- 'THDS'
area <- 'goa'

# query data 
abc <- q_specs(year=year, species=species, area=area, db=db, save=F) |> 
  filter(!area_label == "GOA") |> 
  rename(strata = area_label, ABC = acceptable_biological_catch, Year = year) |> 
  mutate(strata = ifelse(strata == "W", "WGOA", 
                         ifelse(strata == "C", "CGOA", "EGOA")))

DBI::dbDisconnect(db)
 
# abc <- read.csv(paste0(dat_path, "/sst_harvest_specs.csv")) |> 
#   filter(Type == "ABC", Year > 2009) |> 
#   rename(WGOA = W, CGOA = C, EGOA = E) |> 
#   pivot_longer(cols = c(WGOA, CGOA, EGOA), names_to = "strata", values_to = "ABC")
  

db <- DBI::dbConnect (odbc::odbc(),
                      dsn = "afsc",
                      uid = keyring::key_list("afsc")$username,
                      pwd =  keyring::key_get("afsc", keyring::key_list("afsc")$username))


try <- goa_thornyhead(year = 2024)
area_catch <- read.csv(paste0("data/2024_sept_pt/Groundfish Total Catch.csv")) |>
  rename(strata = Regulatory.Area) |> 
  mutate(strata = ifelse(strata == "Western Gulf", "WGOA", 
                         ifelse(strata == "Central Gulf", "CGOA", "EGOA"))) |> 
  group_by(Year, strata) |> 
  summarize(Catch = sum(Catch..mt.))

catch_plot <- left_join(abc, area_catch)

ggplot(catch_plot |> filter(Year > 2009, Year < 2025)) + 
  geom_point(aes(Year, Catch)) +
  geom_line(aes(Year, Catch)) +
  geom_point(aes(Year, ABC), col = "red3") +
  geom_line(aes(Year, ABC), col = "red3") +
  facet_wrap(~factor(strata, levels = c("WGOA", "CGOA", "EGOA")))  +
  ylab("Catch (t)") +
  theme_bw()

ggsave(filename = paste0(out_path, '/catch_by_area.png'),
       dpi = 300, bg = 'white', units = 'in', height = 4, width = 10)

fish_catch <- read.csv(paste0(dat_path, "/Groundfish Total Catch by Fishery.csv")) |>
  group_by(Year, Gear) |> 
  summarize(Catch = sum(Catch..mt.))

ggplot(fish_catch) +
  geom_col(aes(Year, Catch, fill = Gear)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_fishery.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)

target_catch <- read.csv(paste0(dat_path, "/Groundfish Total Catch by Fishery.csv")) |>
  rename(Target = Trip.Target.Group) |> 
  group_by(Year, Target) |> 
  summarize(Catch = sum(Catch..mt.))

ggplot(target_catch) +
  geom_col(aes(Year, Catch, fill = Target)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/catch_by_target.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)



retention_catch <- read.csv(paste0(dat_path, "/Groundfish Total Catch by Fishery.csv")) |>
  rename(Retention = Retained.Discarded) |> 
  group_by(Year, Retention) |> 
  summarize(Catch = sum(Catch..mt.))
  

ggplot(retention_catch) +
  geom_col(aes(Year, Catch, fill = Retention)) +
  ylab("Catch (mt)")

ggsave(filename = paste0(out_path, '/fishery_retention.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5.5, width = 6)
