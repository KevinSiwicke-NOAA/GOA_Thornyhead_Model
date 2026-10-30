# GOA SST biomass estimation using the bottom trawl and longline survey indices

# Model naming conventions:
# Model 24.2a = m24.2: 1 process error 1 q, xtra observation error for BTS and LLS, but the 'a' is because rema changed how extra obs error was calculated
# Model 26, combine 0-500m and 501-700m to accommodate trawl survey restratification. Drops 1990, 93, 96, and 2001 from trawl data


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
dat_path <- paste0("data/", YEAR, "_nov_pt"); dir.create(dat_path)
out_path <- paste0("results/", YEAR, "_nov_pt"); dir.create(out_path)

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
  write_csv(paste0(out_path, '/biomass_2strat_data_wide.csv'))

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
input <- prepare_rema_input(model_name = 'M24.2a',
                            multi_survey = 1,
                            biomass_dat = biomass_dat_3 |> filter(year < 2025),
                            cpue_dat = cpue_dat,
                            sum_cpue_index = TRUE,
                            start_year = 1990,
                            end_year = 2027,
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


# Model 26 (M24.2f) -----
# 26 (formally 24.2f), combination of 0-500 and 501-700 that drops years when 
input <- prepare_rema_input(model_name = 'M26',
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

m26 <- fit_rema(input)
out26 <- tidy_rema(m26)
out26$parameter_estimates

xtra <- tidy_extra_cv(out26)
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

ggsave(filename = paste0(out_path, '/M26_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 7, width = 8)

# compare with addition of 2025
compare <- compare_rema_models(rema_models = list(m24.2a, m26))

compare$plots$biomass_by_cpue_strata  +
  facet_wrap(~factor(strata, levels=c('WGOA', 'CGOA', 'EGOA')), ncol = 3) +
  geom_line(linewidth = 1.2) +
  # geom_point(aes(y = obs), col = 'black') +
  # geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci), col = 'black') +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#D55E00', '#0072B2')) + 
  scale_color_discrete(type = c('#D55E00', '#0072B2')) + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 55000)) +
  theme(legend.position = "top")

ggsave(filename = paste0(out_path, '/M24.2a_26_strata_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 6, width = 10)

compare$plots$total_predicted_biomass +
  theme(legend.position = 'top', legend.text = element_text(size = 12)) +
  geom_line(linewidth = 1.2) +
  labs(x = NULL, y = 'Biomass (t)',
       fill = NULL, colour = NULL, shape = NULL, lty = NULL) +
  scale_fill_discrete(type = c('#D55E00', '#0072B2')) + 
  scale_color_discrete(type = c('#D55E00', '#0072B2')) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 115000))

ggsave(filename = paste0(out_path, '/M24.2a_26_total_biomass.png'),
       dpi = 400, bg = 'white', units = 'in', height = 5, width = 8)

# apportionment ----
# First look at the confidence interval of the regions
# mc_params <- out26$biomass_by_cpue_strata |> filter(year == YEAR + 1) |> 
#   mutate(brd = 0.0225 * pred,
#          sigma = (log(pred_uci) - log(pred_lci)) / (2*1.96))
# 
# log_mus <- mc_params$log_pred
# log_sigmas <- mc_params$sigma
# 
# # Covariance from historical data
# # cov_matrix <- diag(log_sigmas^2)
# wide_log_data <- out26$biomass_by_cpue_strata %>% select(year, strata, log_pred) |> 
#   pivot_wider(names_from = strata, values_from = log_pred) %>%
#   dplyr::select(-year) # Drop the year column to keep only numeric values
# 
# cov_matrix <- cov(wide_log_data, use = "pairwise.complete.obs")
# 
# n_sim = 1000000
# 
# log_samples <- MASS::mvrnorm(n = n_sim, mu = log_mus, Sigma = cov_matrix)
# 
# # 2. Exponentiate to convert to true Lognormal space
# real_samples <- exp(log_samples)
# 
# # 3. Compute row-wise totals for each simulation run
# totals <- rowSums(real_samples)
# 
# # 4. Compute the relative proportion for each group
# proportions <- real_samples / totals
# 
# results <- matrix(NA, nrow = 3, ncol = 3, 
#                   dimnames = list(c(mc_params$strata[1], mc_params$strata[2], mc_params$strata[3]), 
#                                   c("Median_Proportion", "2.5% Lower Bound", "97.5% Upper Bound")))
# for (i in 1:3) {
#   # We extract the 2.5th and 97.5th percentiles of the calculated proportions
#   results[i, ] <- quantile(proportions[, i], probs = c(0.5, 0.025, 0.975))
# }
# 
# # Display your results
# print(round(results, 4))
# 
# for (i in 1:3) {
#   # We extract the 2.5th and 97.5th percentiles of the calculated proportions
#   results[i, ] <- quantile(proportions[, i], probs = c(0.5, 0.025, 0.975))
# }
# 
# # Display your results
# print(round(results, 4))

appo <- compare$output$biomass_by_strata |>
  filter(model_name %in% c("M26")) |> 
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
  write_csv(paste0(out_path, '/M26_apportionment.csv'))

ggplot(appo, aes(year, proportion_std)) + 
  geom_col(aes(fill = strata)) + 
  facet_wrap(~model_name) +
  coord_flip() +
  labs(x = NULL, y = 'Proportion', fill = 'Region')

ggsave(filename = paste0(out_path, '/M26_appo_std.png'),
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
new <- out26$total_predicted_biomass |> filter(year == 2025) |> pull(pred)
(new-old)/old

# Table for predicted biomass by strata with total LCI/UCI
strata_table <- compare$output$biomass_by_strata |>
  filter(model_name %in% c("M26")) |> 
  mutate(strata = ifelse(grepl('CGOA', strata), 'CGOA',
                         ifelse(grepl('EGOA', strata), 'EGOA',
                                'WGOA'))) |>
  group_by(model_name, year, strata) |>
  summarise(biomass = sum(pred, na.rm = T), lci = sum(pred_lci, na.rm = T), uci = sum(pred_uci, na.rm = T)) |>
  pivot_wider(id_cols = c(year, model_name), names_from = strata, values_from = biomass) |> 
  mutate(WGOA = round(WGOA, 0), CGOA = round(CGOA, 0), EGOA = round(EGOA, 0))

write.csv(strata_table, paste0(out_path, '/M26_Table_strata_biomass_pred.csv'))

goa_table <- compare$output$total_predicted_biomass |> 
  filter(model_name %in% c("M26")) |> 
  select(model_name, year, pred, pred_lci, pred_uci) |> 
  mutate(GOA_Total = round(pred, 0), LCI = round(pred_lci, 0), UCI = round(pred_uci, 0)) 

write.csv(goa_table, paste0(out_path, '/M26_Table_goa_biomass_pred.csv'))

# Pull length data and make figures
source("r/Length_figures.r")
source("r/catch_surv_maps.r")

