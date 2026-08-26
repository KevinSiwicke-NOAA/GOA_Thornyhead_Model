# GOA SST biomass estimation using the bottom trawl and longline survey indices

# Look at Model 22 with 3 area process errors compared to a shared single 
# process error.

# Model naming conventions:
# Model 24.1 = m22: accepted in 2022, 3 process errors, 1 q, xtra observation error for BTS and LLS
# Model 24.2 = m24: 1 process error 1 q, xtra observation error for BTS and LLS


# Set up ----

# assessment year
YEAR <- 2026

# Consider whether the rema package needs to be 'updated'
# install.packages("devtools")
# devtools::install_github("JaneSullivan-NOAA/rema", dependencies = TRUE)

libs <- c('rema', 'readr', 'dplyr', 'tidyr', 'ggplot2', 'cowplot')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)

# folder set up
dat_path <- paste0("data/", YEAR, "_sep_pt"); dir.create(dat_path)
out_path <- paste0("results/", YEAR, "_sep_pt"); dir.create(out_path)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# Read data ----
source("r/biom_data_pull.r")

# bottom trawl survey
# biomass_dat <- read_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
# biomass_dat <- read_csv(paste0(dat_path, "/goa_sst_biomass_", YEAR, ".csv"))
biomass_dat_3 <- model_dat$biomass_dat_3strat

biomass_dat_3 |> 
  write_csv(paste0(dat_path, "/goa_sst_biomass_3strat_", YEAR, ".csv"))
biomass_dat_3 |> 
  tidyr::expand(year = min(biomass_dat_3$year):(YEAR),
                strata) |> 
  left_join(biomass_dat_3 |> 
              mutate(value = ifelse(is.na(biomass), NA,
                                    paste0(prettyNum(round(biomass, 0), big.mark = ','), ' (',
                                           format(round(cv, 3), nsmall = 3, trim = TRUE), ')')))) |> 
  pivot_wider(id_cols = c(year), names_from = strata, values_from = value) |> 
  arrange(year) |> 
  write_csv(paste0(out_path, '/biomass_3strat_data_wide.csv'))

biomass_dat_2 <- model_dat$biomass_dat_2strat

biomass_dat_2 |> 
  write_csv(paste0(dat_path, "/goa_sst_biomass_2strat_", YEAR, ".csv"))

biomass_dat_2 |> 
  tidyr::expand(year = min(biomass_dat_2$year):(YEAR),
                strata) |> 
  left_join(biomass_dat_2 |> 
              mutate(value = ifelse(is.na(biomass), NA,
                                    paste0(prettyNum(round(biomass, 0), big.mark = ','), ' (',
                                           format(round(cv, 3), nsmall = 3, trim = TRUE), ')')))) |> 
  pivot_wider(id_cols = c(year), names_from = strata, values_from = value) |> 
  arrange(year) |> 
  write_csv(paste0(out_path, '/biomass_3strat_data_wide.csv'))

# imputed data 
biomass_dat_3_imp <- model_dat$biomass_dat_3strat_impute
biomass_dat_2_imp <- model_dat$biomass_dat_2strat_impute

# longline survey rpws
# cpue_dat <- read_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))
cpue_dat <- model_dat$cpue_dat 
cpue_dat |> 
  write_csv(paste0(dat_path, "/goa_sst_rpw_", YEAR, ".csv"))
cpue_dat |> 
  tidyr::expand(year = min(biomass_dat_3$year):(YEAR),
                strata) |> 
  left_join(cpue_dat |> 
              mutate(value = ifelse(is.na(cpue), NA,
                                    paste0(prettyNum(round(cpue, 0), big.mark = ','), ' (',
                                           format(round(cv, 3), nsmall = 3, trim = TRUE), ')')))) |> 
  pivot_wider(id_cols = c(year), names_from = strata, values_from = value) |> 
  arrange(year) |> 
  write_csv(paste0(out_path, '/cpue_data_wide.csv'))

# Model 24.2a -----
# Same as M24.2 but with rema updates to handling of extra cv (no longer an option because restratification not separated at 500-m depth)
input <- prepare_rema_input(model_name = 'M24.2a_3Base',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_3 |> filter(year < 2025),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2a <- fit_rema(input)
out24.2a <- tidy_rema(m24.2a)
out24.2a$parameter_estimates

xtra <- tidy_extra_cv(out24.2a)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_500','CGOA_0_500','EGOA_0_500',
                                                         'WGOA_501_700','CGOA_501_700','EGOA_501_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                   coord_cartesian(ylim = c(0, 55000)),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "v") 

ggsave(filename = paste0(out_path, '/M24.2a_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Get biomass estimates from M24.2a to compare to the old M24.2
m24.2a_pred <- out24.2a$total_predicted_biomass |> 
  select(model_name, year, pred, pred_lci, pred_uci) |> 
  mutate(GOA_Total = round(pred, 0), LCI = round(pred_lci, 0), UCI = round(pred_uci, 0)) |> 
  filter(year < 2024)

m24.2_pred <- read.csv("results/2024_nov_pt/M24.2_Table_goa_biomass_pred.csv") |> 
  filter(year < 2024) |> 
  mutate(pred_new = m24.2a_pred$pred,
         per_change = 100*(pred_new - pred) / pred)

ggplot(m24.2_pred, aes(year, per_change)) +
  geom_line(linewidth = 2) + 
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = 'Year', y = '% change') + 
  scale_y_continuous(limits = c(-2, 2))

mean(m24.2_pred$per_change)

# Model 24.2b -----
# same as 24.2a but with years missing 501-700, values are imputed (just stepping stone not an option)
input <- prepare_rema_input(model_name = 'M24.2b_3Impute',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_3_imp |> filter(year < 2025),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2b <- fit_rema(input)
out24.2b <- tidy_rema(m24.2b)
out24.2b$parameter_estimates

xtra <- tidy_extra_cv(out24.2b)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_500','CGOA_0_500','EGOA_0_500',
                                                         'WGOA_501_700','CGOA_501_700','EGOA_501_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2b_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Model 24.2c -----
# Same as 24.2a but drops data from years where 501-700 missing samples (just stepping stone not an option)
input <- prepare_rema_input(model_name = 'M24.2c_3Drop',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_3 |> filter(!year %in% c(1990, 1993, 1996, 2001, 2025)),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1992,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 1, 2, 2, 2, 3, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2c <- fit_rema(input)
out24.2c <- tidy_rema(m24.2c)
out24.2c$parameter_estimates

xtra <- tidy_extra_cv(out24.2c)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_500','CGOA_0_500','EGOA_0_500',
                                                         'WGOA_501_700','CGOA_501_700','EGOA_501_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2c_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Compare all models with 9 strata. None can be used, but to get a sense of change from base before including 2025.
compare <- compare_rema_models(rema_models = list(m24.2b, m24.2c, m24.2a))

cowplot::plot_grid(ggplot(data = compare$output$biomass_by_strata, aes(x = year)) +
                     theme(legend.position = 'top', legend.text = element_text(size = 16)) +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_500','CGOA_0_500','EGOA_0_500',
                                                         'WGOA_501_700','CGOA_501_700','EGOA_501_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci, fill = model_name), col = NA, alpha = 0.2) +
                     geom_line(aes(y = pred, col = model_name), linewidth = 1.2) +
                     geom_point(aes(y = obs, fill = model_name), pch = 21) +
                     geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci, col = model_name)) +
                     labs(x = NULL, y = 'Biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) + 
                     scale_color_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) +
                     scale_y_continuous(expand = c(0,0)) +
                     coord_cartesian(ylim = c(0, 55000)), # Only constrains the visual window
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3) +
                     geom_line(linewidth = 1.2) +
                     geom_point(aes(x = year, y = obs), col = 'black') +
                     geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), col = 'black') +
                     labs(x = NULL, y = 'Relative Population Weights',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) + 
                     scale_color_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) + 
                     scale_y_continuous(expand = c(0,0)) +
                     coord_cartesian(ylim = c(0, 42000)) + # Only constrains the visual window
                     theme(legend.position = "none"),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "v")

ggsave(filename = paste0(out_path, '/M24.2abc_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 10, width = 12)

compare$plots$total_predicted_biomass +
  theme(legend.position = 'top', legend.text = element_text(size = 12)) +
  geom_line(linewidth = 1.2) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) + 
  scale_color_discrete(type = c('#009E73', '#E69F00', '#56B4E9')) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 110000))

ggsave(filename = paste0(out_path, '/M24.2abc_total_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5, width = 8)

# Look at the contribution of going imputation and drop relative to base
comp_3 <- bind_rows(out24.2a$total_predicted_biomass, out24.2b$total_predicted_biomass,
                    out24.2c$total_predicted_biomass) |> 
  pivot_wider(id_cols = c(year), names_from = model_name, values_from = pred) |> 
  mutate(M24.2b_Impute = 100 * (M24.2b_3Impute - M24.2a_3Base) / M24.2a_3Base,
         M24.2c_Drop = 100 * (M24.2c_3Drop - M24.2a_3Base) / M24.2a_3Base) |> 
  select(year, M24.2b_Impute, M24.2c_Drop) |> 
  filter (year < 2024) |> 
  pivot_longer(cols = c(M24.2b_Impute, M24.2c_Drop), names_to = 'Model', values_to = 'Per_change')

ggplot(comp_3, aes(year, Per_change, col = Model)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = 'Year', y = '% Change') +
  geom_line(linewidth = 2) +
  scale_color_discrete(type = c('#E69F00', '#56B4E9')) + 
  ggtitle('From M24.2a (Base), 3-strata') +
  theme(legend.position = 'top')

ggsave(filename = paste0(out_path, '/3_strat_change.png'),
       dpi = 400, bg = 'white', units = 'in', height = 4, width = 6)

# Get values for last ten years for doc (2014 to 2023)
comp_3 |> filter(year > 2013) |> group_by(Model) |> summarize(mn_per = mean(Per_change))

# Model 24.2d -----
# same as24.2a but with 1-500 and 501-700 combined naively with no fill in
input <- prepare_rema_input(model_name = 'M24.2d_2Naive',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2 |> filter(year < 2025),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2d <- fit_rema(input)
out24.2d <- tidy_rema(m24.2d)
out24.2d$parameter_estimates

xtra <- tidy_extra_cv(out24.2d)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     coord_cartesian(ylim = c(0, 70000)),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2d_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Model 24.2e -----
# Combination of 0-500 and 501-700 that includes imputed data from 24.2b 
input <- prepare_rema_input(model_name = 'M24.2e_2Impute',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2_imp |> filter(year < 2025),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2e <- fit_rema(input)
out24.2e <- tidy_rema(m24.2e)
out24.2e$parameter_estimates

xtra <- tidy_extra_cv(out24.2e)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     coord_cartesian(ylim = c(0, 70000)),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "v") 

ggsave(filename = paste0(out_path, '/M24.2e_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Model 24.2f -----
# Combination of 0-500 and 501-700 that drops years like 24.2c
input <- prepare_rema_input(model_name = 'M24.2f_2Drop',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2 |> filter(!year %in% c(1990, 1993, 1996, 2001, 2025)),
                            cpue_dat = cpue_dat |> filter(year < 2025),
                            sum_cpue_index = TRUE,
                            start_year = 1992,
                            end_year = 2025,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2f <- fit_rema(input)
out24.2f <- tidy_rema(m24.2f)
out24.2f$parameter_estimates

xtra <- tidy_extra_cv(out24.2f)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2f_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# no data change, just strata change
compare <- compare_rema_models(rema_models = list(m24.2d, m24.2e, m24.2f))

cowplot::plot_grid(ggplot(data = compare$output$biomass_by_strata, aes(x = year)) +
  theme(legend.position = 'top', legend.text = element_text(size = 16)) +
  facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                      'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci, fill = model_name), col = NA, alpha = 0.2) +
  geom_line(aes(y = pred, col = model_name), linewidth = 1.2) +
  geom_point(aes(y = obs, fill = model_name), pch = 21) +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci, col = model_name)) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(ylim = c(0, 71000)) + # Only constrains the visual window
  scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')),
  compare$plots$cpue_by_strata  +
    facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3) +
    geom_line(linewidth = 1.2) +
    geom_point(aes(y = obs), col = 'black') +
    geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), col = 'black') +
    labs(x = NULL, y = 'Relative Population Weights',
         fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
    scale_fill_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) + 
    scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) + 
    theme(legend.position = "none"),
  ncol = 1,
  rel_heights = c(2, 1),
  align = "v")

ggsave(filename = paste0(out_path, '/M24.2def_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 10, width = 12)

compare <- compare_rema_models(rema_models = list(m24.2a, m24.2d, m24.2e, m24.2f))

compare$plots$total_predicted_biomass +
  theme(legend.position = 'top', legend.text = element_text(size = 12)) +
  geom_line(linewidth = 1.2) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#009E73', '#F0E442', '#D55E00', '#0072B2')) + 
  scale_color_discrete(type = c('#009E73', '#F0E442', '#D55E00', '#0072B2')) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 115000))

ggsave(filename = paste0(out_path, '/M24.2adef_total_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5, width = 8)

# Look at the contribution of going from 3 to 2 strata
comp_3_to_2 <- bind_rows(out24.2a$total_predicted_biomass, out24.2b$total_predicted_biomass,
                         out24.2c$total_predicted_biomass, out24.2d$total_predicted_biomass,
                         out24.2e$total_predicted_biomass, out24.2f$total_predicted_biomass) |> 
  pivot_wider(id_cols = c(year), names_from = model_name, values_from = pred) |> 
  mutate(M24.2d_Naive = 100 * (M24.2d_2Naive - M24.2a_3Base) / M24.2a_3Base,
         M24.2e_Impute = 100 * (M24.2e_2Impute - M24.2b_3Impute) / M24.2b_3Impute,
         M24.2f_Drop = 100 * (M24.2f_2Drop - M24.2c_3Drop) / M24.2c_3Drop) |> 
  select(year, M24.2d_Naive, M24.2e_Impute, M24.2f_Drop) |> 
  pivot_longer(cols = c(M24.2d_Naive, M24.2e_Impute, M24.2f_Drop), names_to = 'Model', values_to = 'Per_change')

ggplot(comp_3_to_2 |> filter (year < 2024), aes(year, Per_change, col = Model)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = 'Year', y = '% Change') +
  geom_line(linewidth = 2) +
  scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) +
  ggtitle('From 3 to 2 trawl depth strata') +
  theme(legend.position = 'top')

ggsave(filename = paste0(out_path, '/2_strat_change.png'),
       dpi = 400, bg = 'white', units = 'in', height = 4.5, width = 6)

comp_3_to_2 |> filter(year > 2013) |> group_by(Model) |> summarize(mean_per_change = mean(Per_change)) # Both are ~ <1% increase overall

comp_2_to_base <- bind_rows(out24.2a$total_predicted_biomass,
                         out24.2e$total_predicted_biomass, out24.2f$total_predicted_biomass) |> 
  pivot_wider(id_cols = c(year), names_from = model_name, values_from = pred) |> 
  mutate(M24.2e_Impute = 100 * (M24.2e_2Impute - M24.2a_3Base) / M24.2a_3Base,
         M24.2f_Drop = 100 * (M24.2f_2Drop - M24.2a_3Base) / M24.2a_3Base) |> 
  select(year, M24.2e_Impute, M24.2f_Drop) |> 
  pivot_longer(cols = c(M24.2e_Impute, M24.2f_Drop), names_to = 'Model', values_to = 'Per_change')

comp_2_to_base |> filter(year > 2013) |> group_by(Model) |> summarize(mean_per_change = mean(Per_change)) # Both are ~ <1% increase overall

# Updated with 2025 data
# Model 24.2d_25 -----
# same as 24.2d with 1-500 and 501-700 combined naively with no fill in and added 2025 data
input <- prepare_rema_input(model_name = 'M24.2d_25_2Naive',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2027,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2d_25 <- fit_rema(input)
out24.2d_25 <- tidy_rema(m24.2d_25)
out24.2d_25$parameter_estimates

xtra <- tidy_extra_cv(out24.2d_25)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     coord_cartesian(ylim = c(0, 70000)),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2d_25_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Model 24.2e_25 -----
# 24.2e updated with 2025 data
input <- prepare_rema_input(model_name = 'M24.2e_25_2Impute',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2_imp,
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2027,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2e_25 <- fit_rema(input)
out24.2e_25 <- tidy_rema(m24.2e_25)
out24.2e_25$parameter_estimates

xtra <- tidy_extra_cv(out24.2e_25)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     coord_cartesian(ylim = c(0, 70000)),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2e_25_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# Model 24.2f_25 -----
# 24.2f, combination of 0-500 and 501-700 that drops years when 
input <- prepare_rema_input(model_name = 'M24.2f_25_2Drop',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_2 |> filter(!year %in% c(1990, 1993, 1996, 2001)),
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1992,
                            end_year = 2027,
                            PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
                            q_options = list(
                              pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
                              pointer_q_cpue = c(1, 1, 1)),
                            extra_biomass_cv = list(assumption = 'extra_cv'),
                            extra_cpue_cv = list(assumption = 'extra_cv'))

m24.2f_25 <- fit_rema(input)
out24.2f_25 <- tidy_rema(m24.2f_25)
out24.2f_25$parameter_estimates

xtra <- tidy_extra_cv(out24.2f_25)
xtra$biomass_by_strata <- xtra$biomass_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

xtra$cpue_by_strata <- xtra$cpue_by_strata |> 
  mutate(tot_cv = sqrt(exp(tot_sd_log_obs^2) - 1),
         delta_cv = tot_cv - obs_cv)

cvplots <- plot_extra_cv(xtra)
cowplot::plot_grid(cvplots$biomass_by_strata +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3),
                   cvplots$cpue_by_strata +
                     labs(y = 'RPW') +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "h") 

ggsave(filename = paste0(out_path, '/M24.2f_25_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# compare with addition of 2025
compare <- compare_rema_models(rema_models = list(m24.2d_25, m24.2e_25, m24.2f_25))

cowplot::plot_grid(ggplot(data = compare$output$biomass_by_strata, aes(x = year)) +
                     theme(legend.position = 'top', legend.text = element_text(size = 16)) +
                     facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
                                                         'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
                     geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci, fill = model_name), col = NA, alpha = 0.2) +
                     geom_line(aes(y = pred, col = model_name), linewidth = 1.2) +
                     geom_point(aes(y = obs, fill = model_name), pch = 21) +
                     geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci, col = model_name)) +
                     labs(x = NULL, y = 'Biomass (t)',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) +
                     scale_y_continuous(expand = c(0,0)) +
                     coord_cartesian(ylim = c(0, 71000)) + # Only constrains the visual window
                     scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')),
                   compare$plots$cpue_by_strata  +
                     facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3) +
                     geom_line(linewidth = 1.2) +
                     geom_point(aes(y = obs), col = 'black') +
                     geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), col = 'black') +
                     labs(x = NULL, y = 'Relative Population Weights',
                          fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
                     scale_fill_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) + 
                     scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) + 
                     theme(legend.position = "none"),
                   ncol = 1,
                   rel_heights = c(2, 1),
                   align = "v")

ggsave(filename = paste0(out_path, '/M24.2def_25_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 10, width = 12)

compare <- compare_rema_models(rema_models = list(m24.2a, m24.2d_25, m24.2e_25, m24.2f_25))

compare$plots$total_predicted_biomass +
  theme(legend.position = 'top', legend.text = element_text(size = 12)) +
  geom_line(linewidth = 1.2) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#009E73','#F0E442', '#D55E00', '#0072B2')) + 
  scale_color_discrete(type = c('#009E73', '#F0E442', '#D55E00', '#0072B2')) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 115000))

ggsave(filename = paste0(out_path, '/M24.2adef_25_total_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5, width = 8)

# And compare that to the contribution from adding 2025
comp_new_data <- bind_rows(out24.2d$total_predicted_biomass, out24.2d_25$total_predicted_biomass,
                           out24.2e$total_predicted_biomass, out24.2e_25$total_predicted_biomass,
                         out24.2f$total_predicted_biomass, out24.2f_25$total_predicted_biomass) |> 
  pivot_wider(id_cols = c(year), names_from = model_name, values_from = pred) |> 
  mutate(M24.2d_Naive = 100 * (M24.2d_25_2Naive - M24.2d_2Naive) / M24.2d_2Naive,
         M24.2e_Impute = 100 * (M24.2e_25_2Impute - M24.2e_2Impute) / M24.2e_2Impute,
         M24.2f_Drop = 100 * (M24.2f_25_2Drop - M24.2f_2Drop) / M24.2f_2Drop) |> 
  select(year, M24.2d_Naive, M24.2e_Impute, M24.2f_Drop)  |> 
  pivot_longer(cols = c(M24.2d_Naive, M24.2e_Impute, M24.2f_Drop), names_to = 'Option', values_to = 'Per_change')
  
ggplot(comp_new_data, aes(year, Per_change, col = Option)) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = 'Year', y = '% Change') +
  geom_line(linewidth = 2) +
  scale_color_discrete(type = c('#F0E442', '#D55E00', '#0072B2')) +
  ggtitle("From '23 to '25, 2 trawl depth strata") +
  theme(legend.position = 'top')

ggsave(filename = paste0(out_path, '/2025_change.png'),
       dpi = 400, bg = 'white', units = 'in', height = 4.5, width = 6)

# Not really for 2026, but look at putting a prior or PE, and not estimating
# two additional observation errors as they lead to model validation issues
# Model 24.2g_25 -----
# 24.2g, same as f, but only additional cv on longline, none on biomass
# input <- prepare_rema_input(model_name = 'M24.2g_25_prior_xtra',
#                             multi_survey = 1,
#                             biomass_dat = biomass_dat_2 |> filter(!year %in% c(1990, 1993, 1996, 2001)),
#                             cpue_dat = cpue_dat,
#                             sum_cpue_index = TRUE,
#                             start_year = 1992,
#                             end_year = 2027,
#                             PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1),
#                                                    penalty_options = 'normal_prior',
#                                                    penalty_values = c(-2.3, 0.5)),
#                             q_options = list(
#                               pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
#                               pointer_q_cpue = c(1, 1, 1)),
#                             extra_cpue_cv = list(assumption = 'extra_cv'))
# 
# m24.2g_25 <- fit_rema(input)
# out24.2g_25 <- tidy_rema(m24.2g_25)
# out24.2g_25$parameter_estimates
# 
# input <- prepare_rema_input(model_name = 'M24.2h_25_prior',
#                             multi_survey = 1,
#                             biomass_dat = biomass_dat_2 |> filter(!year %in% c(1990, 1993, 1996, 2001)),
#                             cpue_dat = cpue_dat,
#                             sum_cpue_index = TRUE,
#                             start_year = 1992,
#                             end_year = 2027,
#                             PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1),
#                                               penalty_options = 'normal_prior',
#                                               penalty_values = c(-2.3, 0.1)),
#                             q_options = list(
#                               pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
#                               pointer_q_cpue = c(1, 1, 1)))
# 
# m24.2h_25 <- fit_rema(input)
# out24.2h_25 <- tidy_rema(m24.2h_25)
# out24.2h_25$parameter_estimates
# 
# input <- prepare_rema_input(model_name = 'M24.2i_25_cv',
#                             multi_survey = 1,
#                             biomass_dat = biomass_dat_2 |> filter(!year %in% c(1990, 1993, 1996, 2001)) |> 
#                               mutate(cv = ifelse(cv < 0.15, 0.15, cv)),
#                             cpue_dat = cpue_dat |> 
#                               mutate(cv = ifelse(cv < 0.2, 0.2, cv)),
#                             sum_cpue_index = TRUE,
#                             start_year = 1992,
#                             end_year = 2027,
#                             PE_options = list(pointer_PE_biomass = c(1, 1, 1, 1, 1, 1)),
#                             q_options = list(
#                               pointer_biomass_cpue_strata = c(1, 1, 2, 2, 3, 3),
#                               pointer_q_cpue = c(1, 1, 1)))
# 
# m24.2i_25 <- fit_rema(input)
# out24.2i_25 <- tidy_rema(m24.2i_25)
# out24.2i_25$parameter_estimates
# 
# # 
# # # compare with addition of 2025
# compare <- compare_rema_models(rema_models = list(m24.2f_25, m24.2h_25, m24.2i_25))
# 
# cowplot::plot_grid(ggplot(data = compare$output$biomass_by_strata, aes(x = year)) +
#                      theme(legend.position = 'top', legend.text = element_text(size = 16)) +
#                      facet_wrap(~factor(strata, levels=c('WGOA_0_700','CGOA_0_700','EGOA_0_700',
#                                                          'WGOA_701_1000','CGOA_701_1000','EGOA_701_1000')), ncol = 3) +
#                      geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci, fill = model_name), col = NA, alpha = 0.2) +
#                      geom_line(aes(y = pred, col = model_name), linewidth = 1.2) +
#                      geom_point(aes(y = obs, fill = model_name), pch = 21) +
#                      geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci, col = model_name)) +
#                      labs(x = NULL, y = 'Biomass (t)',
#                           fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
#                      scale_fill_discrete(type = c('#0072B2', 'black', 'firebrick')) +
#                      scale_y_continuous(expand = c(0,0)) +
#                      coord_cartesian(ylim = c(0, 71000)) + # Only constrains the visual window
#                      scale_color_discrete(type = c('#0072B2', 'black', 'firebrick')),
#                    compare$plots$cpue_by_strata  +
#                      facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3) +
#                      geom_line(linewidth = 1.2) +
#                      geom_point(aes(y = obs), col = 'black') +
#                      geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci, col = model_name)) +
#                      labs(x = NULL, y = 'Relative Population Weights',
#                           fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
#                      scale_fill_discrete(type = c('#0072B2', 'black', 'firebrick')) +
#                      scale_color_discrete(type = c('#0072B2', 'black', 'firebrick')) +
#                      theme(legend.position = "none"),
#                    ncol = 1,
#                    rel_heights = c(2, 1),
#                    align = "v")
# 
# # compare <- compare_rema_models(rema_models = list(m24.2a, m24.2d_25, m24.2e_25, m24.2f_25))
# 
# compare$plots$total_predicted_biomass +
#   theme(legend.position = 'top', legend.text = element_text(size = 12)) +
#   geom_line(linewidth = 1.2) +
#   labs(x = NULL, y = 'Biomass (t)',
#        fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
#   scale_fill_discrete(type = c('#0072B2', 'black', 'firebrick')) + 
#   scale_color_discrete(type = c('#0072B2', 'black', 'firebrick')) +
#   scale_y_continuous(expand = c(0,0), limits = c(0, 115000))

# apportionment ----
appo <- compare$output$biomass_by_strata |>
  filter(model_name %in% c("M24.2d_25_2Naive", "M24.2e_25_2Impute", "M24.2f_25_2Drop")) |> 
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) |>
  group_by(model_name, strata, year) |>
  summarize(stratum_biomass = sum(pred)) |>
  group_by(model_name, year) |>
  mutate(total_biomass = sum(stratum_biomass)) |>
  ungroup() |>
  mutate(proportion_std = stratum_biomass / total_biomass) |>
  mutate(strata = factor(strata, labels = c('EGOA', 'CGOA', 'WGOA'), levels = c('EGOA', 'CGOA', 'WGOA'), ordered = TRUE)) |>
  arrange(year, strata)

appo |>
  pivot_wider(id_cols = c(strata, year), names_from = model_name, values_from = proportion_std) |>
  write_csv(paste0(out_path, '/M24.2d_apportionment.csv'))

ggplot(appo, aes(year, proportion_std)) + 
  geom_col(aes(fill = strata)) + 
  facet_wrap(~model_name) +
  coord_flip() +
  labs(x = NULL, y = 'Proportion', fill = 'Region')

ggsave(filename = paste0(out_path, '/M24.2f_appo_std.png'),
       dpi = 600, bg = 'white', units = 'in', height = 5, width = 10)

full_sumtable <- appo |>
  filter(year == 2027) |>
  mutate(natmat = 0.03,
         OFL = natmat * total_biomass,
         maxABC = 0.75 * natmat * total_biomass,
         ABC = maxABC)

sumtable <- full_sumtable |>
  distinct(model_name, year, biomass = total_biomass, OFL, maxABC) |>
  select(model_name, year, biomass, OFL, maxABC) |>
  write_csv(paste0(out_path, '/abc_ofl_summary.csv'))

# percent changes -----
biomass_dat_2 |> filter(year %in% c(2023,2025)) |>
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) |>
  group_by(year, strata) |>
  summarise(biomass = sum(biomass, na.rm = T)) |>
  pivot_wider(id_cols = c(strata), names_from = year, values_from = biomass) |>
  mutate(percent_change = (`2025`-`2023`)/`2023`)

cpue_dat |> filter(year %in% c(2023, 2025)) |>
  pivot_wider(id_cols = c(strata), names_from = year, values_from = cpue) |>
  mutate(percent_change = (`2025`-`2023`)/`2023`)

old <- 59459
new <- out24.2f_25$total_predicted_biomass |> filter(year == 2025) |> pull(pred)
(new-old)/old

# Table for predicted biomass by strata with total LCI/UCI
strata_table <- compare$output$biomass_by_strata |>
  filter(model_name %in% c("M24.2d_25_2Naive", "M24.2e_25_2Impute", "M24.2f_25_2Drop")) |> 
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) |>
  group_by(model_name, year, strata) |>
  summarise(biomass = sum(pred, na.rm = T), lci = sum(pred, na.rm = T), uci = sum(pred, na.rm = T)) |>
  pivot_wider(id_cols = c(year, model_name), names_from = strata, values_from = biomass) |> 
  mutate(WGOA = round(WGOA, 0), CGOA = round(CGOA, 0), EGOA = round(EGOA, 0))

write.csv(strata_table, paste0(out_path, '/M24.2def_Table_strata_biomass_pred.csv'))

goa_table <- compare$output$total_predicted_biomass |> 
  filter(model_name %in% c("M24.2d_25_2Naive", "M24.2e_25_2Impute", "M24.2f_25_2Drop")) |> 
  select(model_name, year, pred, pred_lci, pred_uci) |> 
  mutate(GOA_Total = round(pred, 0), LCI = round(pred_lci, 0), UCI = round(pred_uci, 0)) 

write.csv(goa_table, paste0(out_path, '/M24.2def_Table_goa_biomass_pred.csv'))

# Pull length data and make figures
source("r/Length_figures.r")

