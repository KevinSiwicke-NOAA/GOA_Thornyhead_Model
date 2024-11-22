#### GOA THORNYHEAD SARA FILE ####
# Directory Setup ----
YEAR <- 2024

# folder set up
dat_path <- paste0("data/2024_nov_pt")
out_path <- paste0("results/2024_nov_pt")

SARAdir <- paste0(out_path, "/SARA")

if(file.exists(SARAdir)==F) {
  dir.create(paste0(SARAdir), showWarnings = T)
}

## Files and data sets used in SARA file
harvest_specs <- read.csv( paste0(results_path, "rema_biom_OROX_SC.csv"))
hist_harvest <- read.csv(paste0(dir_main, "Data\\static\\OROX_historical_harvestlimits.csv"))

# rema estimated biomass
# note that the rema biomass for combined T4/5 is different from T4 + T5 
## Use combined 4/5 if I need CIs, otherwise use separate Tier 4 & 5 (as in assessment) 
##  to get total B 
#T45_rema_biom <- read.csv( paste0(results_path, "rema_biom_OROX_T45_new.csv")) %>%
#  filter( area == "GOA" )
T4_rema_biom <- read.csv( paste0(results_path, "rema_biom_OROX_SC.csv")) %>%
  filter( area == "GOA" )
T5_rema_biom <- read.csv( paste0(results_path, "rema_biom_OROX_T5_new.csv")) %>%
  filter( area == "GOA" )
rema_biom <- full_join(T4_rema_biom, T5_rema_biom) %>%
  group_by(year) %>%
  summarise(biomass = round(sum(biomass, na.rm=T ), 0 ) )

rema_biom_AYR <- rema_biom %>%
  filter(year == max(year)) %>%
  mutate(lci = "NA", uci = "NA") # add this in case decided to include lci/ uci numbers


# trawl biomass for Tier 4 and 5
surv_biom <- read.csv(paste0(raw_data_path, "race_goa_biomass_raw_", AYR, ".csv" ) ) %>%
  left_join(orox_rf %>% select(species_name, species_code= race_code, tier= new_tier)) %>%
  filter(year >= 1990, tier %in% c(4,5)) %>%
  group_by(year) %>%
  summarise(surv_biom = round(sum(area_biomass, na.rm=T)))

# what cindy called the race data OROXbiomass_dat
# catch data
catch_sara_dat <- read.csv(paste0(dat_path, "/Groundfish Total Catch.csv")) %>%
  filter(year >= 2003 ) %>%
  filter(!(cas_code == 136 & fmp_subarea %in% c("CG", "WG"))) %>% # removing northerns
  filter(!(cas_code %in% c(cas_codes_dsr) & fmp_subarea == "SE")) %>% # removing dsr from SE
  group_by(year) %>%
  summarise(catch = round(sum(tot_catch, na.rm=T))) %>%
  arrange(year)


#### Static SARA info ----
stock <- "SST"
stock_name <- "GULF OF ALASKA THORNYHEAD ROCKFISH COMPLEX"
fmp <- "GOA"
AYR <- 2024 #assessment year
asmt_type <- "Operational full"
asmt_mo <- "Dec"
tier <- "5"
num_sexes <- "NA"
num_fish <- "NA"
rec_mult <- "NA"
recage <- "NA"
complex <- "Yes"
asmt_mod_cat <- "2- Index-Based; NA"
asmt_mod <- "Index Method"
mod_version <- "24.2"
lead_lab <- "AFSC"
email <- "kevin.siwicke@noaa.gov"
review_result <- "Full acceptance" #
catch_input_dat <- 0 #options: 0 - None, 1 - Major gaps preclude use, 2 - Major gaps in some sectors(s), 3 - Minor gaps across sectors, 4 - Minor gaps in some sector(s), 5 - Near complete knowledge)
abund_input_dat <- 3 # options: 0 - None, 1 - Uncertain or expert opinion, 2 - Standardized fishery-dependent, 3 - Limited fishery-independent, 4 - Comprehensive fishery-independent, 5 - Absolute abundance)
biol_input_dat <- 2 # options: 0 - None, 1 - Proxy-based, 2 - Empirical and proxy-based, 3 - Mostly empirical estimates, 4 - Track changes over time, 5 - Comprehensive over time and space)
sizeage_comp_input_dat <- 0 # options: 0 - None, 1 - Major gaps preclude use, 2 - Support data-limited only, 3 - Gaps, bus supports age-structured assessment, 4 - Support fishery composition, 5 - Very complete)
ecosys_link <- 0 

# F section
Fbasis <- "Catch/Total Stock Biomass"
Funit <- "Metric tons"
Flim_basis <- "F35% as proxy; F = M; max catch"  # trying to figure out if this is FOFL or FABC
Fmsy_basis <- "F35% as proxy; F = M; max catch" # trying to figure out if this is FOFL or FABC

# Biomass section
b_basis <- "Total Tier 4 & 5 Biomass" 
b_unit <- "Metric Tons"
est_method <- "NA" # (options Tiers 1-3 only: Asymptotic, Credible, Bootstrapped, "NA" for Tiers 4-6)
interval_size <- "NA" # Specify size of confidence interval (options Tiers 1-3 only: 50 to 99, "NA" for Tiers 4-6)
b_msy_basis <- "NA" # (options Tiers 1-3 only: Direct estimate, S_MSY escapement, Average Survey CPUE, B40%, B35%, B30%, "NA" for Tiers 4-6)

# numbers
Fest <- catch_sara_dat %>% filter(year %in% c(AYR-1)) %>% pull(catch)  # last full year's catch for OROX (do last exploitation for others)
Flim <- hist_harvest %>% filter(Year == AYR & HL == "OFL") %>% pull(Total) # OFL spec year (not new rec spec)
Fmsy <- hist_harvest %>% filter(Year == AYR & HL == "OFL") %>% pull(Total) # OFL last year

b_msy <- "NA" # (Tiers 1-3 only, "NA" for Tiers 4-6)
surv_desc <- "GOA bottom trawl survey" # ("NA" for Tier 6)

stock_status <- "SAFE report in 2023 indicates that this stock was not subjected to overfishing or is being overfished."

# Other fishery numbers that are not used in this file
ages <- "NA"
recruit <- "NA"
spawn_biom <- "NA"
F_yrs <- "NA"


# Compile the .dat file ----
cat(
  "#STOCK","\n",
  stock, "\n",
  "#STOCK_NAME","\n",
  stock_name, "\n",
  "#REGION","\n",
  fmp, "\n",
  "#ASMT_TYPE", "\n",
  asmt_type, "\n",
  "#ASMT_YEAR", "\n",
  AYR, "\n",
  "#ASMT_MONTH", "\n",
  asmt_mo, "\n",
  "#TIER", "\n",
  tier, "\n",
  "#NUM_SEXES", "\n",
  num_sexes, "\n",
  "#NUM_FISHERIES", "\n",
  num_fish, "\n",
  "#REC_MULT", "\n",
  rec_mult, "\n",
  "#REC_MULT", "\n",
  recage, "\n",
  "#COMPLEX", "\n",
  complex, "\n",
  "#LAST_DATA_YEAR", "\n",
  rema_biom_AYR$year, "\n",
  "ASMT_MODEL_CATEGORY", "\n",
  asmt_mod_cat, "\n",
  "#ASMT_MODEL", "\n",
  asmt_mod, "\n", 
  "#MODEL_VERSION", "\n",
  mod_version, "\n",
  "#ENSEMBLE", "\n",
  "NA", "\n",
  "#LEAD_LAB", "\n",
  lead_lab, "\n",
  "#POC_EMAIL", "\n",
  email, "\n",
  "#REVIEW_RESULT", "\n",
  review_result, "\n",
  "#CATCH_INPUT_DATA", "\n",
  catch_input_dat, "\n",
  "#ABUNDANCE_INPUT_DATA", "\n",
  abund_input_dat, "\n",
  "#BIOLOGICAL_INPUT_DATA", "\n",
  biol_input_dat, "\n",
  "#SIZEAGE_COMP_INPUT_DATA", "\n",
  sizeage_comp_input_dat, "\n",
  "#ECOSYSTEM_LINKAGE", "\n",
  ecosys_link, 
  "#FISHING_MORTALITY_ESTIMATES ----------------------------------------------------------------------------------------", "\n",
  "F_YEAR", "\n",
  AYR-1, "\n",
  "#F_BASIS", "\n",
  Fbasis, "\n",
  "#F_UNIT", "\n",
  Funit, "\n",
  "#BEST_F_ESTIMATE", "\n",
  Fest, "\n",
  "#F_LIMIT", "\n",
  Flim, "\n",
  "#F_LIMIT_BASIS", "\n",
  Flim_basis, "\n",
  "#F_MSY", "\n",
  Fmsy, "\n",
  "#F_MSY_BASIS", "\n",
  Fmsy_basis,
  "#BIOMASS_ESTIMATES --------------------------------------------------------------------------------------------------", "\n",
  "#B_YEAR", "\n",
  rema_biom_AYR$year, "\n",
  "#B_BASIS", "\n",
  b_basis, "\n",
  "#B_UNIT", "\n",
  b_unit, "\n",
  "#BEST_B_ESTIMATE", "\n",
  rema_biom_AYR$biomass, "\n",
  "#LOWER_B_ESTIMATE", "\n",
  rema_biom_AYR$lci, "\n",
  "#UPPER_B_ESTIMATE", "\n",
  rema_biom_AYR$uci, "\n",
  "#ESTIMATE_METHOD", "\n",
  est_method, "\n",
  "#INTERVAL_SIZE", "\n",
  interval_size, "\n",
  "#B_MSY", "\n",
  b_msy, "\n",
  "#B_MSY_BASIS", "\n",
  b_msy_basis,
  "#TIME_SERIES_ESTIMATES  ----------------------------------------------------------------------------------------------", "\n",
  "#FISHERYYEAR", "\n",
  rema_biom$year, "\n",
  "#AGE", "\n",
  ages, "\n",
  "#RECRUITMENT", "\n",
  recruit, "\n",
  "#SPAWNBIOMASS", "\n",
  spawn_biom, "\n",
  "#TOTALBIOMASS", "\n",
  rema_biom$biomass, "\n",
  "#TOTFSHRYMORT", "\n",
  F_yrs, "\n",
  "#TOTALCATCH", "\n",
  catch_sara_dat$catch, "\n",
  "#SURVEYDESC", "\n",
  surv_desc, "\n",
  "#STOCKNOTES", "\n",
  stock_status,
  "#SURVEY_ESTIMATES [OPTIONAL] ------------------------------------------------------------------------------------------", "\n",
  surv_biom$surv_biom, "\n", 
  file=paste0(SARAdir,"/OROCKGOA",AYR,".dat"))
  
# End