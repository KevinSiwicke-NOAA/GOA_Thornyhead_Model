# GOA SST biomass estimation using the bottom trawl and longline survey indices

# Look at Model 22 with 3 area process errors compared to a shared single 
# process error.

# Model naming conventions:
# Model 24.1 = m22: accepted in 2022, 3 process errors, 1 q, xtra observation error for BTS and LLS
# Model 24.2 = m24: 1 process error 1 q, xtra observation error for BTS and LLS


# Set up ----

# assessment year
YEAR <- 2023

# Consider whether the rema package needs to be 'updated'
# install.packages("devtools")
# devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)

libs <- c('rema', 'readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# folder set up
dat_path <- paste0("data/2024_nov_pt"); dir.create(dat_path)
out_path <- paste0("results/2024_nov_pt"); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# Read data ----
source("r/biom_data_pull.r")

# bottom trawl survey
# biomass_dat <- read_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
# biomass_dat <- read_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
biomass_dat <- model_dat$biomass_dat |> 
  filter(year < (YEAR + 1))

biomass_dat |> 
  write_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
biomass_dat |> 
  tidyr::expand(year = min(biomass_dat$year):(YEAR),
                strata) |> 
  left_join(biomass_dat |> 
              mutate(value = ifelse(is.na(biomass), NA,
                                    paste0(prettyNum(round(biomass, 0), big.mark = ','), ' (',
                                           format(round(cv, 3), nsmall = 3, trim = TRUE), ')')))) |> 
  pivot_wider(id_cols = c(year), names_from = strata, values_from = value) |> 
  arrange(year) |> 
  write_csv(paste0(out_path, '/biomass_data_wide.csv'))

# longline survey rpws
# cpue_dat <- read_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))
cpue_dat <- model_dat$cpue_dat |> 
  filter(year < (YEAR + 1))
cpue_dat |> 
  write_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))
  
# Model 24.1_3_PE -----
# Same as M22, but updated log transformation
input <- prepare_rema_input(model_name = 'Model 24.1_3pe',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = YEAR + 2,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 2, 2, 2, 3, 3, 3)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.1 <- fit_rema(input)
out24.1 <- tidy_rema(m24.1)
out24.1$parameter_estimates

# Model 24.2_1_PE -----
# Same as M24.1 but single PE
input <- prepare_rema_input(model_name = 'Model 24.2_1pe',
                            multi_survey = 1,
                            biomass_dat = biomass_dat,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = YEAR + 2,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2 <- fit_rema(input)
out24.2 <- tidy_rema(m24.2)
out24.2$parameter_estimates
# Compare estimated extra obs err (error) for 3-PE and 1-PE ----
compare <- compare_rema_models(rema_models = list(m24.1, m24.2))

cowplot::plot_grid(compare$plots$biomass_by_strata +
                     theme(legend.position = 'top', legend.text = element_text(size = 16)) +
                     facet_wrap(~factor(strata, levels=c('WGOA (0-500 m)','WGOA (501-700 m)','WGOA (701-1000 m)',
                                                         'CGOA (0-500 m)','CGOA (501-700 m)','CGOA (701-1000 m)',
                                                         'EGOA (0-500 m)','EGOA (501-700 m)','EGOA (701-1000 m)')), ncol = 3) +
                     geom_line(linewidth = 1.2) +
                     labs(x = NULL, y = 'Biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#440154FF', '#2A788EFF')) +
                     scale_color_discrete(type = c('#440154FF', '#2A788EFF')),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 1) +
                     geom_line(linewidth = 1.2) +
                     labs(x = NULL, y = 'Relative Population Weights',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#440154FF', '#2A788EFF')) +
                     scale_color_discrete(type = c('#440154FF', '#2A788EFF')) +
                     theme(legend.position = "none"),
                   ncol = 2,
                   rel_widths = c(2.5, 1),
                   align = "h")

ggsave(filename = paste0(out_path, '/M24.1_3pe_vs_M24.2_1pe_fits.png'),
       dpi = 400, bg = 'white', units = 'in', height = 9, width = 14)

compare$plots$total_predicted_biomass +
  theme(legend.position = 'top', legend.text = element_text(size = 16)) +
  geom_line(linewidth = 1.2) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#440154FF', '#2A788EFF')) +
  scale_color_discrete(type = c('#440154FF', '#2A788EFF'))

ggsave(filename = paste0(out_path, '/M24.1_3pe_vs_M24.2_1pe_totalbiomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 3.5, width = 8)

# param estimates
compare$output$parameter_estimates %>%
  write_csv(paste0(out_path, '/M24.1_3pe_M24.2_1pe_parameters.csv'))

# predicted biomass by strata and total for each model
biom <- compare$output$biomass_by_strata %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred) %>%
  bind_rows(compare$output$total_predicted_biomass %>%
              mutate(strata = 'Total') %>%
              pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = pred)) %>%
  write_csv(paste0(out_path, '/M24.1_3pe_M24.2_1pe_biomass_pred.csv'))

# apportionment ----
appo <- compare$output$biomass_by_strata %>%
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) %>%
  group_by(model_name, strata, year) %>%
  summarize(stratum_biomass = sum(pred)) %>%
  group_by(model_name, year) %>%
  mutate(total_biomass = sum(stratum_biomass)) %>%
  ungroup() %>%
  mutate(proportion_std = stratum_biomass / total_biomass) %>%
  mutate(strata = factor(strata, labels = c('EGOA', 'CGOA', 'WGOA'), levels = c('EGOA', 'CGOA', 'WGOA'), ordered = TRUE)) %>%
  arrange(year, strata)

appo %>%
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = proportion_std) %>%
  write_csv(paste0(out_path, '/M24.1_3pe_M24.2_1pe_apportionment.csv'))

ggplot(appo, aes(year, proportion_std)) + 
  geom_col(aes(fill = strata)) + 
  facet_wrap(~model_name) +
  coord_flip() +
  labs(x = NULL, y = 'Proportion', fill = 'Region')

ggsave(filename = paste0(out_path, '/m24.1_m24.2_appo_std.png'),
       dpi = 600, bg = 'white', units = 'in', height = 5, width = 10)

full_sumtable <- appo %>%
  filter(year > (YEAR-2)) %>%
  mutate(natmat = 0.03,
         OFL = natmat * total_biomass,
         maxABC = 0.75 * natmat * total_biomass,
         ABC = maxABC)

sumtable <- full_sumtable %>%
  distinct(model_name, year, biomass = total_biomass, OFL, maxABC) %>%
  select(model_name, year, biomass, OFL, maxABC) %>%
  write_csv(paste0(out_path, '/abc_ofl_summary.csv'))

# percent changes -----
biomass_dat %>% filter(year %in% c(2021,2023)) %>%
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) %>%
  group_by(year, strata) %>%
  summarise(biomass = sum(biomass, na.rm = T)) %>%
  pivot_wider(id_cols = c(strata), names_from = year, values_from = biomass) %>%
  mutate(percent_change = (`2023`-`2021`)/`2021`)

cpue_dat %>% filter(year %in% c(2023, 2022)) %>%
  pivot_wider(id_cols = c(strata), names_from = year, values_from = cpue) %>%
  mutate(percent_change = (`2023`-`2022`)/`2022`)

old <- out24.1$total_predicted_biomass %>% filter(year == 2024) %>% pull(pred)
new <- out24.2$total_predicted_biomass %>% filter(year == 2024) %>% pull(pred)
(new-old)/old


# Table for predicted biomass by strata with total LCI/UCI
strata_table <- compare$output$biomass_by_strata %>%
  filter(model_name == "Model 24.2_1pe") %>% 
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) %>%
  group_by(year, strata) %>%
  summarise(biomass = sum(pred, na.rm = T), lci = sum(pred, na.rm = T), uci = sum(pred, na.rm = T) ) %>%
  pivot_wider(id_cols = c(year), names_from = strata, values_from = biomass) %>% 
  mutate(WGOA = round(WGOA, 0), CGOA = round(CGOA, 0), EGOA = round(EGOA, 0)) %>% 
  write.csv(paste0(out_path, '/M24.2_1pe_Table_strata_biomass_pred.csv'))

goa_table <- compare$output$total_predicted_biomass %>% 
  filter(model_name == "Model 24_1pe") %>% 
  select(year, pred, pred_lci, pred_uci) %>% 
  mutate(GOA_Total = round(pred, 0), LCI = round(pred_lci, 0), UCI = round(pred_uci, 0)) %>% 
  write.csv(paste0(out_path, '/M24_1pe_Table_goa_biomass_pred.csv'))

# Pull length data and make figures
source("r/Length_figures.r")
