library(tmbstan)

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
# run simulation testing
par_ests_24.2f <- sim_test(mod_name = m24.2f_25, replicates = 500, cpue = TRUE) # note some warnings
length(which(is.na(par_ests_24.2f[,1])))
mod_par_ests_f <- data.frame("log_PE" = m24.2f_25$env$last.par[1],
                           "log_q" = m24.2f_25$env$last.par[2],
                           "log_tau_biomass" = m24.2f_25$env$last.par[3],
                           "log_tau_cpue" = m24.2f_25$env$last.par[4],
                           type = "model")
names(par_ests_24.2f) <- names(mod_par_ests_f) # rename for ease
thrn_par_ests_f <- rbind(mod_par_ests_f, par_ests_24.2f) # recovered and model parameters in one data frame
thrn_par_ests_f$sp <- rep("GOA Thornyhead")

thrn_sim_f <- thrn_par_ests_f %>% pivot_longer(1:4, names_to = "parameter")

# Relative Error: ((om-em)/om)
sim_re_f <- thrn_sim_f %>%
  pivot_wider(id_cols = c("sp", "parameter"),
              names_from = type, values_from = value) %>%
  unnest(cols = c(model, recovered)) %>%
  mutate(RE = (model-recovered)/model*100) %>%
  group_by(sp, parameter) %>%
  mutate(label = paste0("Median RE=", formatC(median(RE, na.rm = TRUE), format = "f", digits = 1), "%")) %>%
  suppressWarnings()

p1thrn <- plot_sim(thrn_sim_f, "GOA Thornyhead Simulation")

p2thrn <- plot_re(sim_re_f %>% filter(sp == "GOA Thornyhead"))

cowplot::plot_grid(p1thrn, p2thrn, ncol = 1)
ggsave(filename = paste0(out_path, '/sim_f.png'), width = 11, height = 7, units = "in", bg = "white")

par_ests_24.2h <- sim_test(mod_name = m24.2h_25, replicates = 500, cpue = TRUE)
length(which(is.na(par_ests_24.2h[,1])))
mod_par_ests_h <- data.frame("log_PE" = m24.2h_25$env$last.par[1],
                             "log_q" = m24.2h_25$env$last.par[2],
                             type = "model")
names(par_ests_24.2h) <- names(mod_par_ests_h) # rename for ease
thrn_par_ests_h <- rbind(mod_par_ests_h, par_ests_24.2h) # recovered and model parameters in one data frame
thrn_par_ests_h$sp <- rep("GOA Thornyhead")

thrn_sim_h <- thrn_par_ests_h %>% pivot_longer(1:2, names_to = "parameter")

# Relative Error: ((om-em)/om)
sim_re_h <- thrn_sim_h %>%
  pivot_wider(id_cols = c("sp", "parameter"),
              names_from = type, values_from = value) %>%
  unnest(cols = c(model, recovered)) %>%
  mutate(RE = (model-recovered)/model*100) %>%
  group_by(sp, parameter) %>%
  mutate(label = paste0("Median RE=", formatC(median(RE, na.rm = TRUE), format = "f", digits = 1), "%")) %>%
  suppressWarnings()

p1thrn <- plot_sim(thrn_sim_h, "GOA Thornyhead Simulation")
p2thrn <- plot_re(sim_re_h %>% filter(sp == "GOA Thornyhead"))
cowplot::plot_grid(p1thrn, p2thrn, ncol = 1)
ggsave(filename = paste0(out_path, '/sim_h.png'), width = 11, height = 7, units = "in", bg = "white")

par_ests_24.2i <- sim_test(mod_name = m24.2i_25, replicates = 500, cpue = TRUE)
length(which(is.na(par_ests_24.2i[,1])))
mod_par_ests_i <- data.frame("log_PE" = m24.2i_25$env$last.par[1],
                             "log_q" = m24.2i_25$env$last.par[2],
                             type = "model")
names(par_ests_24.2i) <- names(mod_par_ests_i) # rename for ease
thrn_par_ests_i <- rbind(mod_par_ests_i, par_ests_24.2i) # recovered and model parameters in one data frame
thrn_par_ests_i$sp <- rep("GOA Thornyhead")

thrn_sim_i <- thrn_par_ests_i %>% pivot_longer(1:2, names_to = "parameter")

# Relative Error: ((om-em)/om)
sim_re_i <- thrn_sim_i %>%
  pivot_wider(id_cols = c("sp", "parameter"),
              names_from = type, values_from = value) %>%
  unnest(cols = c(model, recovered)) %>%
  mutate(RE = (model-recovered)/model*100) %>%
  group_by(sp, parameter) %>%
  mutate(label = paste0("Median RE=", formatC(median(RE, na.rm = TRUE), format = "f", digits = 1), "%")) %>%
  suppressWarnings()

p1thrn <- plot_sim(thrn_sim_i, "GOA Thornyhead Simulation")
p2thrn <- plot_re(sim_re_i %>% filter(sp == "GOA Thornyhead"))
cowplot::plot_grid(p1thrn, p2thrn, ncol = 1)
ggsave(filename = paste0(out_path, '/sim_i.png'), width = 11, height = 7, units = "in", bg = "white")


# thrn_osa <- rema::get_osa_residuals(thrn_mod)
# thrn_osa$plots$qq
# thrn_osa$plots$biomass_qq
# thrn_osa$plots$cpue_qq
# thrn_osa$plots$biomass_resids
# thrn_osa$plots$cpue_resids

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
library(tmbstan)
thrn_comp_log_tau_f <- mcmc_comp(mod_name = m24.2f_25, it_numb = 5000, chain_numb = 3) 
thrn_comp_log_tau_h <- mcmc_comp(mod_name = m24.2h_25, it_numb = 5000, chain_numb = 3) 
thrn_comp_log_tau_i <- mcmc_comp(mod_name = m24.2i_25, it_numb = 5000, chain_numb = 3) 

# clean up data frame names by renaming things
thrn_qq <- thrn_comp_log_tau_f[[1]]
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
m <- summary(thrn_comp_log_tau_f[[3]])$summary[1:4, c(6,4,8,10), drop = FALSE]
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
posterior_draws <- as.array(thrn_comp_log_tau_f[[3]])
library(bayesplot)
mcmc_pairs(posterior_draws, pars = c("log_PE", "log_q", "log_tau_biomass", "log_tau_cpue"),
           diag_fun = "hist",    # Puts histograms on the diagonal
           off_diag_fun = "scatter")  # Puts scatterplots on the off-diagonal

rstan::traceplot(thrn_comp_log_tau_f[[3]], ncol = 3) + ggtitle(label = NULL, subtitle = "GOA Thornyhead")
ggsave(paste0("vignettes/ex4_traceplots.png"), width = 11, height = 8, units = "in", bg = "white")



# clean up data frame names by renaming things
thrn_qq <- thrn_comp_log_tau_h[[1]]
thrn_qq$par_name <- recode(thrn_qq$par,
                           V1 = "log_PE",
                           V2 = "log_q")
thrn_qq$sp <- rep("GOA Thornyhead")
qq_dat <- thrn_qq

plot_laplace_mcmc(thrn_qq, "GOA Thornyhead")
# ggsave(paste0("vignettes/ex4_mcmc_qq.png"), width = 11, height = 7, units = "in", bg = "white")

# parameter summary table - compare the marginal max likelihood estimates (MMLE)
# with the mcmc estimates. Rhat is the potential scale reduction factor on split
# chains (at convergence, Rhat=1)
m <- summary(thrn_comp_log_tau_f[[3]])$summary[1:4, c(6,4,8,10), drop = FALSE]
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
posterior_draws <- as.array(thrn_comp_log_tau_h[[3]])
mcmc_pairs(posterior_draws, pars = c("log_PE", "log_q"),
           diag_fun = "hist",    # Puts histograms on the diagonal
           off_diag_fun = "scatter")  # Puts scatterplots on the off-diagonal

rstan::traceplot(thrn_comp_log_tau_h[[3]], ncol = 3) + ggtitle(label = NULL, subtitle = "GOA Thornyhead")
# ggsave(paste0("vignettes/ex4_traceplots.png"), width = 11, height = 8, units = "in", bg = "white")

# clean up data frame names by renaming things
thrn_qq <- thrn_comp_log_tau_i[[1]]
thrn_qq$par_name <- recode(thrn_qq$par,
                           V1 = "log_PE",
                           V2 = "log_q" )
thrn_qq$sp <- rep("GOA Thornyhead")
qq_dat <- thrn_qq

plot_laplace_mcmc(thrn_qq, "GOA Thornyhead")
# ggsave(paste0("vignettes/ex4_mcmc_qq.png"), width = 11, height = 7, units = "in", bg = "white")

# parameter summary table - compare the marginal max likelihood estimates (MMLE)
# with the mcmc estimates. Rhat is the potential scale reduction factor on split
# chains (at convergence, Rhat=1)
m <- summary(thrn_comp_log_tau_f[[3]])$summary[1:4, c(6,4,8,10), drop = FALSE]
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
posterior_draws <- as.array(thrn_comp_log_tau_i[[3]])
library(bayesplot)
mcmc_pairs(posterior_draws, pars = c("log_PE", "log_q"),
           diag_fun = "hist",    # Puts histograms on the diagonal
           off_diag_fun = "scatter")  # Puts scatterplots on the off-diagonal

rstan::traceplot(thrn_comp_log_tau_i[[3]], ncol = 3) + ggtitle(label = NULL, subtitle = "GOA Thornyhead")
ggsave(paste0("vignettes/ex4_traceplots.png"), width = 11, height = 8, units = "in", bg = "white")