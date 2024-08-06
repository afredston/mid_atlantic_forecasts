fit_drm <- function(amarel = FALSE,
                    use_poisson_link = FALSE,
                    run_name = "test",
                    results_path = file.path("results", run_name),
                    create_dir = TRUE,
                    do_dirichlet = 1,
                    eval_l_comps = 0,
                    T_dep_mortality = 0,
                    T_dep_recruitment = 0,
                    T_dep_movement = 0,
                    spawner_recruit_relationship = 1,
                    run_forecast = 0,
                    process_error_toggle = 1,
                    exp_yn = 0,
                    warmup = 1000,
                    iter = 2000,
                    max_treedepth =  10,
                    chains =  1,
                    refresh = 10,
                    cores = 1,
                    adapt_delta = 0.95, 
                    drm_name = "process_sdm",
                    number_quantiles = 3,
                    quantiles_calc = c(0.05, 0.5, 0.95),
                    known_f = 0,
                    known_historic_f = 1,
                    #          cv = 0.1, # replaced by pr_sigma_obs_sigma, below 
                    h = 0.8,
                    dcap = 1 / 3, 
                    # all priors that Stan needs for parameter estimation 
                    pr_init_dep_alpha = 1.5,
                    pr_init_dep_beta = 3,
                    pr_beta_obs_mu = 0.001,
                    pr_beta_obs_sigma = 0.1,
                    pr_beta_obs_int_mu = -100,
                    pr_beta_obs_int_sigma = 4,
                    pr_raw_mu = 0,
                    pr_raw_sigma = 1,
                    pr_sigma_r_raw_mu = .2,
                    pr_sigma_r_raw_sigma = .1,
                    pr_sigma_obs_mu = 0.2, 
                    pr_sigma_obs_sigma = 0.2, # formerly "cv" / "sigma_obs_cv"
                    pr_log_mean_recruits_mu = 7,
                    pr_log_mean_recruits_sigma = 5,
                    pr_log_r0_mu = 15,
                    pr_log_r0_sigma = 5,
                    pr_Topt_mu = 18 ,
                    pr_Topt_sigma = 2,
                    pr_width_mu = 4,
                    pr_width_sigma = 2,
                    pr_beta_t_mu = 0,
                    pr_beta_t_sigma = 2,
                    pr_beta_rec_mu = 0,
                    pr_beta_rec_sigma = 2,
                    pr_alpha_alpha = 12, # concentrated around 0.4
                    pr_alpha_beta = 20, # concentrated around 0.4
                    pr_d_mu = 0.1, # diffusion rate as a proportion of total population size within the patch
                    pr_d_sigma = 0.1,
                    pr_p_length_50_sel_sigma = .2, # note that there is no mu prior for this parameter passed in because it's calculated in the code from length_50_sel_guess and loo 
                    pr_sel_delta_mu = 2,
                    pr_sel_delta_sigma = 4,
                    pr_theta_d_mu = 0.5,
                    pr_theta_d_sigma = 0.5,
                    patch_r0s = 0
) {
  if (amarel == TRUE) {
    dyn.load('/projects/community/gcc/9.2.0/gc563/lib64/libgfortran.so.5')
    
  }
  
  if (create_dir == TRUE) {
    if (!dir.exists(results_path)) {
      dir.create(results_path, recursive = TRUE)
    }
  }
  
  if (amarel == FALSE) {
    load(here::here("processed-data", "stan_data_prep.Rdata"))
  }
  
  if (amarel == TRUE) {
    load(
      '/home/fredston/mid_atlantic_forecasts/processed-data/stan_data_prep.Rdata'
    )
  }
  
  area <- rep(1,np)
  # by default f is a time-varying instantaneous rate 
  # if we want to turn that off: 
  if (known_historic_f == 0) {
    f <- matrix(0.2, nrow = nrow(f), ncol = ncol(f))
    # estimated f in 2017 from NOAA SAW 66
  }
  
  # keep the values from the final year of f for the projections
  if (known_f == 0) {
    f_proj <-
      matrix(rep(f[, ncol(f)], ncol(f_proj)),
             ncol = ncol(f_proj),
             nrow = nrow(f_proj))
  }
  
  if (T_dep_movement == 1) {
    dcap = 1
  }
  stan_data <- list(
    np = np,
    patches = patches,
    n_ages = n_ages,
    ny_train = ny,
    ny_proj = ny_proj,
    n_lbins = n_lbins,
    n_at_length = len,
    dens = dens,
   #  area = area, 
  #  area =  matrix(1, nrow = nrow(dens), ncol = ncol(dens)),
    sbt = sbt,
    sbt_proj = sbt_proj,
    m = m,
    f = f,
    f_proj = f_proj,
    k = k,
    loo = loo,
    t0 = t0,
    #    cv = cv,
    length_50_sel_guess = length_50_sel_guess,
    age_sel = age_sel,
    bin_mids = bin_mids,
    sel_100 = sel_100,
    age_at_maturity = age_at_maturity,
    l_at_a_key = l_at_a_mat,
    wt_at_age = wt_at_age,
    do_dirichlet = do_dirichlet,
    eval_l_comps = eval_l_comps,
    T_dep_mortality = T_dep_mortality,
    T_dep_recruitment = T_dep_recruitment,
    T_dep_movement = T_dep_movement,
    spawner_recruit_relationship = spawner_recruit_relationship,
    run_forecast = run_forecast,
    exp_yn = exp_yn,
    process_error_toggle = process_error_toggle,
    number_quantiles = number_quantiles,
    quantiles_calc = quantiles_calc,
    h = h,
    dcap = dcap,
    use_poisson_link = use_poisson_link,
    pr_init_dep_alpha = pr_init_dep_alpha,
    pr_init_dep_beta = pr_init_dep_beta,
    pr_beta_obs_mu = pr_beta_obs_mu,
    pr_beta_obs_sigma = pr_beta_obs_sigma,
    pr_beta_obs_int_mu = pr_beta_obs_int_mu,
    pr_beta_obs_int_sigma = pr_beta_obs_int_sigma,
    pr_raw_mu = pr_raw_mu,
    pr_raw_sigma = pr_raw_sigma,
    pr_sigma_r_raw_mu = pr_sigma_r_raw_mu,
    pr_sigma_r_raw_sigma = pr_sigma_r_raw_sigma,
    pr_sigma_obs_mu = pr_sigma_obs_mu, 
    pr_sigma_obs_sigma = pr_sigma_obs_sigma,
    pr_log_mean_recruits_mu = pr_log_mean_recruits_mu,
    pr_log_mean_recruits_sigma = pr_log_mean_recruits_sigma,
    pr_log_r0_mu = pr_log_r0_mu,
    pr_log_r0_sigma = pr_log_r0_sigma,
    pr_Topt_mu = pr_Topt_mu ,
    pr_Topt_sigma = pr_Topt_sigma,
    pr_width_mu = pr_width_mu,
    pr_width_sigma = pr_width_sigma,
    pr_beta_t_mu = pr_beta_t_mu,
    pr_beta_t_sigma =pr_beta_t_sigma,
    pr_beta_rec_mu = pr_beta_rec_mu,
    pr_beta_rec_sigma = pr_beta_rec_sigma,
    pr_alpha_alpha = pr_alpha_alpha,
    pr_alpha_beta = pr_alpha_beta,
    pr_d_mu =pr_d_mu,
    pr_d_sigma = pr_d_sigma,
    pr_p_length_50_sel_sigma = pr_p_length_50_sel_sigma,
    pr_sel_delta_mu = pr_sel_delta_mu,
    pr_sel_delta_sigma = pr_sel_delta_sigma,
    pr_theta_d_mu = pr_theta_d_mu,
    pr_theta_d_sigma = pr_theta_d_sigma,
  area = area,
  patch_r0s = patch_r0s
  )
  nums <- 100 * exp(-.2 * (0:(n_ages - 1)))
  
  if (amarel == FALSE) {
    stan_file = here::here("src", paste0(drm_name, ".stan"))
  }
  
  if (amarel == TRUE) {
    stan_file = paste0('/home/fredston/mid_atlantic_forecasts/src/',
                       drm_name,
                       '.stan')
  }
  drm_model <- cmdstan_model( here::here("src",paste0(drm_name, ".stan")))
  
  stan_model_fit = drm_model$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = warmup,
    iter_sampling = iter - warmup,
    parallel_chains = cores,
    refresh = refresh,
    max_treedepth = max_treedepth,
    adapt_delta = 0.85,
    init = lapply(1:chains, function(x)
      list(
        sigma_obs = 0.1,
        Topt = jitter(12, 4),
        log_r0 = jitter(10, 5),
        beta_obs = jitter(1e-6, 4),
        beta_obs_int = jitter(-10, 2)
      ))
  )
  
  stan_model_fit$save_object(file = file.path(results_path,
                                              "stan_model_fit.rds"))
  
  # fit <- readRDS(file.path(results_path,
  #                          "stan_model_fit.rds"))
  # browser()
  #   
  #   stan_model_fit <- stan(
  #     sample_file = file.path(results_path,paste0(run_name,".csv")),
  #     diagnostic_file = file.path(results_path,paste0(run_name,"_diagnostics.csv")),
  #     file = stan_file,
  #     data = stan_data,
  #     chains = chains,
  #     warmup = warmup,
  #     iter = iter,
  #     cores = cores,
  #     refresh = refresh,
  #     control = list(max_treedepth = 8,
  #                    adapt_delta = 0.95),
  #     init = lapply(1:chains, function(x)
  #       list(
  #         Topt = jitter(12, 4),
  #         log_r0 = jitter(10, 5),
  #         beta_obs = jitter(1e-6, 4),
  #         beta_obs_int = jitter(-10, 2)
  #       ))
  #   )
  rm(stan_model_fit)
  # readr::write_rds(stan_model_fit, file = file.path(results_path,
  #                                                 "stan_model_fit.rds"))
  # didn't work on HPC
  
  # readr::write_rds(stan_model_fit, file = file.path(results_path,
  #                                                   "stan_model_fit.rds"))
  #
  gc()
  return("all done")
  
} # close fit_drm