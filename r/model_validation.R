library(tmbstan)

# next, get the thornyhead data set up
thrn_bio_dat <- read.csv(system.file("example_data/goa_thornyhead_2022_biomass_dat.csv", package = "rema"))
thrn_cpue_dat <- read.csv(system.file("example_data/goa_thornyhead_2022_cpue_dat.csv", package = "rema"))

thrn_input <- prepare_rema_input(model_name = 'Model 24.2f_25',
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

# run the model
thrn_mod <- fit_rema(thrn_input) 

# tidy output and plot fitted data
tidy_thrn <- tidy_rema(thrn_mod)
tidy_thrn$parameter_estimates

p2 <- plot_rema(tidy_thrn, biomass_ylab = 'Biomass (t)')$biomass_by_strata + 
  ggtitle(label = "Model Fits to the GOA Thornyhead Data",
          subtitle = "Trawl Survey Biomass Strata") + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = "#21918c", fill = "#21918c", alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
p2
p3 <- plot_rema(tidy_thrn, cpue_ylab = "RPW")$cpue_by_strata + 
  ggtitle(label = NULL, subtitle = "Longline Survey Relative Population Weight (RPW) Strata") + 
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci),
              col = "#21918c", fill = "#21918c", alpha = 0.4) +
  geom_line() +
  geom_point(aes(x = year, y = obs)) +
  geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci))
p3
cowplot::plot_grid(p2, p3, rel_heights = c(0.7, 0.3), ncol = 1)
# ggsave(paste0("vignettes/ex4_thrn_fits.png"), width = 11, height = 10, units = "in", bg = "white")

sim_test <- function(mod_name, replicates, cpue) {
  
  # storage things
  re_est <- matrix(NA, replicates, length(mod_name$par)) # parameter estimates
  
  # go through the model
  suppressMessages(for(i in 1:replicates) {
    
    sim <- mod_name$simulate(complete = TRUE) # simulates the data
    
    # simulated biomass observations:
    tmp_biomass <- matrix(data = exp(sim$log_biomass_obs), ncol = ncol(mod_name$input$data$biomass_obs))
    colnames(tmp_biomass) <- colnames(mod_name$data$biomass_obs)
    # simulated cpue observations, when applicable:
    if (cpue) {tmp_cpue <- matrix(data = sim$cpue_obs, ncol = ncol(mod_name$input$data$cpue_obs))
    colnames(tmp_cpue) <- colnames(mod_name$data$cpue_obs)}
    # set up new data for input
    newinput <- mod_name$input
    newinput$data$biomass_obs <- tmp_biomass # biomass data
    if (cpue) {newinput$data$cpue_obs <- tmp_cpue} # cpue data
    
    # create "obsvec" which is used internally in cpp file as the observation
    # vector for all observation (log biomass + cpue) in the likelihood
    # functions (and required for OSA residuals). note the transpose t() needed
    # to get these matrices in the correct order (by row instead of by col) --
    # note obsvec is masked from users normally in prepare_rema_input()
    newinput$data$obsvec <- t(sim$log_biomass_obs)[!is.na(t(sim$log_biomass_obs))]
    if (cpue) {newinput$data$obsvec <- c(t(sim$log_biomass_obs)[!is.na(t(sim$log_biomass_obs))], t(sim$log_cpue_obs)[!is.na(t(sim$log_cpue_obs))])}
    
    # refit model
    mod_new <- fit_rema(newinput, do.sdrep = FALSE)
    
    # add parameter estimates to matrix
    if(mod_new$opt$convergence == 0) {
      re_est[i, ] <- mod_new$env$last.par[1:length(mod_name$par)]
    } else {
      re_est[i, ] <- rep(NA, length(mod_name$par))
    }
    
  })
  
  re_est <- as.data.frame(re_est); re_est$type <- rep("recovered")
  
  return(re_est)
  
}

# run for thorny and prep data frame -- same process
# run simulation testing
par_ests <- sim_test(mod_name = thrn_mod, replicates = 500, cpue = TRUE) # note some warnings
# sometimes spits out: In stats::nlminb(model$par, model$fn, model$gr, control =
# list(iter.max = 1000,...: NA/NaN function evaluation
n_not_converged_thrn <- length(which(is.na(par_ests[,1]))); n_thrn <- length(par_ests[,1])
prop_converged_thrn <- 1-n_not_converged_thrn/n_thrn
nrow(par_ests)
# get data frame with simulation values for each parameter
mod_par_ests <- data.frame("log_PE" = thrn_mod$env$last.par[1],
                           "log_q" = thrn_mod$env$last.par[2],
                           "log_tau_biomass" = thrn_mod$env$last.par[3],
                           "log_tau_cpue" = thrn_mod$env$last.par[4],
                           type = "model")
names(par_ests) <- names(mod_par_ests) # rename for ease
thrn_par_ests <- rbind(mod_par_ests, par_ests) # recovered and model parameters in one data frame
thrn_par_ests$sp <- rep("GOA Thornyhead")

# reogranize data
thrn_sim <- thrn_par_ests %>% pivot_longer(1:4, names_to = "parameter")
sim_dat <- thrn_sim

# Relative Error: ((om-em)/om)
sim_re <- sim_dat %>%
  pivot_wider(id_cols = c("sp", "parameter"),
              names_from = type, values_from = value) %>%
  unnest(cols = c(model, recovered)) %>%
  mutate(RE = (model-recovered)/model*100) %>%
  group_by(sp, parameter) %>%
  mutate(label = paste0("Median RE=", formatC(median(RE, na.rm = TRUE), format = "f", digits = 1), "%")) %>%
  suppressWarnings()

# plot distribution of simulated parameter estimates
plot_sim <- function(sim_dat, plot_title, fill_col = "#21918c") {
  ggplot(NULL, aes(parameter, value)) +
    # add distribution of recovered parameters
    geom_violin(data = sim_dat %>% filter(type == "recovered"),
                fill = fill_col, alpha = 0.6, draw_quantiles = 0.5) +
    # add "true values" from original model
    geom_point(data = sim_dat %>% filter(type == "model"),
               size = 2, col = "black") +
    facet_wrap(~ parameter, scales = "free", nrow = 1) +
    labs(x = NULL, y = "Parameter estimate", title = plot_title,
         subtitle = "Distribution of parameters estimates (median=horizontal line, true value=point)") +
    scale_x_discrete(labels = NULL, breaks = NULL)
}
# plot relative error
plot_re <- function(sim_re, fill_col = "#21918c") {
  ggplot(NULL, aes(parameter, RE)) +
    # add distribution of recovered parameters
    geom_boxplot(data = sim_re, alpha = 0.6, fill = fill_col,
                 na.rm = TRUE, outlier.size = 0.8) +
    geom_hline(yintercept = 0) +
    # separate by parameter
    facet_wrap(~ parameter+label, scales = "free", nrow = 1) + #, space = "free") +
    labs(x = NULL, y = "Relative error (%)", subtitle = "Distribution of relative error (RE; i.e., (true-estimated values)/true value*100)") +
    scale_x_discrete(labels = NULL, breaks = NULL)
}

p1thrn <- plot_sim(thrn_sim, "GOA Thornyhead Simulation")

p2thrn <- plot_re(sim_re %>% filter(sp == "GOA Thornyhead"))

cowplot::plot_grid(p1thrn, p2thrn, ncol = 1)
# ggsave(paste0("vignettes/ex4_sim_thrn.png"), width = 11, height = 7, units = "in", bg = "white")

thrn_osa <- rema::get_osa_residuals(thrn_mod)
thrn_osa$plots$qq
thrn_osa$plots$biomass_qq
thrn_osa$plots$cpue_qq
thrn_osa$plots$biomass_resids
thrn_osa$plots$cpue_resids

# function to (1) run models with and without laplace approximation for
# comparison and (2) return posterior data frames and models
# function input: model (e.g., AI Pcod or GOA Thornyhead), number of iterations
# (samples), number of chains
mcmc_comp <- function(mod_name, it_numb, chain_numb) {
  
  # set up MCMC chain information
  it_num <- it_numb
  chain_num <- chain_numb
  # mod_name = pcod_mod; it_num = 1000; chain_num = 4
  
  # run model with laplace approximation
  mod_la <- tmbstan(obj = mod_name, chains = chain_num, init = mod_name$par, laplace = TRUE, iter = it_num)
  # run model without laplace approximation, i.e., all parameters fully estimated without assumptions of normality
  mod_mcmc <- tmbstan(obj = mod_name, chains = chain_num, init = mod_name$par, laplace = FALSE, iter = it_num)
  
  # posteriors as data frame
  post_la <- as.data.frame(mod_la); post_la$type <- ("la")
  post_mcmc <- as.data.frame(mod_mcmc); post_mcmc <- post_mcmc[, c(1:length(mod_name$par), dim(post_mcmc)[2])]; post_mcmc$type <- rep("mcmc")
  # informational things... this is for getting the posterior draws, i.e., to
  # test traceplots and such
  post_la$chain <- rep(1:chain_num, each = it_num/2)
  post_la$iter_num <- rep(1:(it_num/2), chain_num)
  post_mcmc$chain <- rep(1:chain_num, each = it_num/2)
  post_mcmc$iter_num <- rep(1:(it_num/2), chain_num)
  post_draws <- rbind(post_la, post_mcmc)
  
  # get quantiles
  qv <- seq(from = 0, to = 1, by = 0.01)
  quant_dat <- data.frame(quant = NULL,
                          la = NULL,
                          mcmc = NULL,
                          par = NULL)
  
  for (i in 1:(dim(post_la)[2]-4)) { # post_la has type, chains, lp, iteration number column that don"t count
    
    tmp <- data.frame(quant = qv,
                      la = quantile(unlist(post_la[i]), probs = qv),
                      mcmc = quantile(unlist(post_mcmc[i]), probs = qv),
                      par = rep(paste0("V", i)))
    
    quant_dat <- rbind(quant_dat, tmp)
    
  }
  
  return(list(quant_dat, post_draws, mod_la, mod_mcmc))
  
}

# run models
thrn_comp_log_tau <- mcmc_comp(mod_name = thrn_mod, it_numb = 5000, chain_numb = 3) 

# clean up data frame names by renaming things
thrn_qq <- thrn_comp_log_tau[[1]]
thrn_qq$par_name <- recode(thrn_qq$par,
                           V1 = "log_PE",
                           V2 = "log_q",
                           V3 = "log_tau_biomass",
                           V4 = "log_tau_cpue" )
thrn_qq$sp <- rep("GOA Thornyhead")
qq_dat <- thrn_qq

# plot
plot_laplace_mcmc <- function(qq_dat, plot_title) {
  ggplot(qq_dat, aes(mcmc, la)) +
    geom_abline(intercept = 0, slope = 1, col = "lightgray") +
    geom_point() +
    facet_wrap(~par_name, scales = "free", nrow = 2,  dir = "v") +
    labs(x = "MCMC", y = "Laplace approx.", title = plot_title) +
    theme(plot.title = element_text(size = 12),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 11),
          strip.text = element_text(size = 11))
}
plot_laplace_mcmc(thrn_qq, "GOA Thornyhead")
# ggsave(paste0("vignettes/ex4_mcmc_qq.png"), width = 11, height = 7, units = "in", bg = "white")

# parameter summary table - compare the marginal max likelihood estimates (MMLE)
# with the mcmc estimates. Rhat is the potential scale reduction factor on split
# chains (at convergence, Rhat=1)
m <- summary(thrn_comp_log_tau[[3]])$summary[1:4, c(6,4,8,10), drop = FALSE]
psum <- bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'MCMC_Laplace', FixedEffect = row.names(m), m, row.names = NULL))
m <- summary(thrn_comp_log_tau[[4]])$summary[1:4, c(6,4,8,10), drop = FALSE]
psum <- psum %>% bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'MCMC_withoutLaplace', FixedEffect = row.names(m), m, row.names = NULL))
pars <- unlist(as.list(thrn_mod$sdrep, "Est"))
pars <- pars[names(pars)[grepl("log_PE|log_q|log_tau", names(pars))]]
sd <- unlist(as.list(thrn_mod$sdrep, "Std"))
sd <- sd[names(sd)[grepl("log_PE|log_q|log_tau", names(sd))]]
psum <- psum %>% bind_rows(data.frame(Stock = 'GOA Thornyhead', Model = 'Base_MMLE_Laplace', FixedEffect = rownames(summary(thrn_comp_log_tau[[3]])$summary[1:4,,drop=FALSE]),
                                      X50. = pars, sd = sd, row.names = NULL) %>%
                             mutate(X2.5. = X50.-1.96*sd,
                                    X97.5. = X50.+1.96*sd,
                                    Rhat = NA) %>%
                             select(-sd))
# write.csv(psum, 'vignettes/ex4_param_table.csv')
posterior_draws <- as.array(thrn_comp_log_tau[[3]])

mcmc_pairs(posterior_draws, pars = c("log_PE", "log_q", "log_tau_biomass", "log_tau_cpue"),
           diag_fun = "hist",    # Puts histograms on the diagonal
           off_diag_fun = "scatter")  # Puts scatterplots on the off-diagonal

p1_mcmc_thrn <- rstan::traceplot(thrn_comp_log_tau[[3]], ncol = 3) + ggtitle(label = NULL, subtitle = "GOA Thornyhead")
ggsave(paste0("vignettes/ex4_traceplots.png"), width = 11, height = 8, units = "in", bg = "white")