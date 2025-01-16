
set.seed(42)
library(tidyverse)
library(ggridges)
library(tidybayes)
#library(Cairo)
library(here)
library(magrittr)
library(rstan)
library(Matrix)
library(rstanarm)
library(cmdstanr)
library(data.table)
library(bayesplot)
set.seed(424242)
rstan_options(javascript = FALSE, auto_write = TRUE)

# load fit_drm function 
# fit_drm() fits the model and writes out the model object and a plot to the results directory
funs <- list.files("functions")
sapply(funs, function(x)
  source(file.path("functions", x)))

# which range edges should be calculated?
#  
quantiles_calc <- c(0.05, 0.5, 0.95)
quantiles_calc <- c(0.95, 0.95, 0.95)


ctrl_file <- read_csv("control_file.csv") %>%
  filter(id == "v0.64")
  
  # filter(
  #   eval_l_comps == 0,
  #   spawner_recruit_relationship == 1,
  #   process_error_toggle == 1,
  #   known_f == 1,
  #   T_dep_mortality == 1
  # ) 
  # ungroup() |>
  # slice(1)

fit_drms <- TRUE
use_poisson_link <- FALSE
if (use_poisson_link){
  run_name <- "yes-pois"
} else {
  run_name <- "no-pois"
}
make_plots <- TRUE
write_summary <- TRUE
iters <- 2000
warmups <- 1000
chains <- 4
cores <- 4

ctrl_file <- ctrl_file |> 
  slice(1)

load(here("processed-data","stan_data_prep.Rdata"))

for(k in 1:nrow(ctrl_file)){
  i = ctrl_file$id[k]  
  
  results_path <- file.path("results",i)
  
  # turn off if you just want to load already-fitted models and analyze them
  
  if (fit_drms==TRUE){
    drm_fits <-  ctrl_file %>%
      filter(id == i)
    
    drm_fits$fits <- list(tryCatch(fit_drm(
      amarel = FALSE,
      use_poisson_link = use_poisson_link,
      create_dir = TRUE,
      run_name = run_name,
      do_dirichlet = drm_fits$do_dirichlet,
      eval_l_comps = drm_fits$eval_l_comps,
      T_dep_movement = drm_fits$T_dep_movement,
      T_dep_mortality = drm_fits$T_dep_mortality,
      T_dep_recruitment = drm_fits$T_dep_recruitment,
      spawner_recruit_relationship = drm_fits$spawner_recruit_relationship,
      process_error_toggle = drm_fits$process_error_toggle,
      exp_yn = drm_fits$exp_yn,
      known_f = drm_fits$known_f,
      known_historic_f = drm_fits$known_historic_f,
      warmup = warmups,
      iter = iters,
      chains = chains,
      cores = cores,
      adapt_delta = 0.85, 
      run_forecast = 1,
      quantiles_calc = quantiles_calc,
      drm_name = "process_sdm",
      patch_r0s = 1
    )
    ) 
    )# as currently written this just adds a column to drm_fits that says "all done". the column `fits` used to contain the model object itself 
    
    
  } # close fit_drms 
}

diagnostic_fit <- read_rds(here("results",run_name, "stan_model_fit.rds"))

diagnostic_fit$diagnostic_summary()

draws <- diagnostic_fit$draws(format = "df")

posterior_drm <- posterior::as_draws_array(diagnostic_fit)

nuts_diagnostics <- bayesplot::nuts_params(diagnostic_fit)

# mcmc_parcoord(
#   posterior_drm,
#   pars = vars(
#     "sigma_obs",
#     "log_r0"
#   ),
#   np = nuts_diagnostics,
#   transformations = scale
# )

# "sigma_r_raw",
# "width",
# "Topt",
# "alpha",
# "sel_delta",
# "p_length_50_sel",
# "beta_obs",
# "theta_d",
# "beta_t",
# "beta_obs_int",

diagnostics <- diagnostic_fit$sampler_diagnostics(format = "df")

diagnostics

mcmc_hist(diagnostic_fit$draws("d"))

test <- tidybayes::spread_draws(
  diagnostic_fit,
  log_r0,
  sigma_obs
) 

# test <- tidybayes::spread_draws(
#   diagnostic_fit,
#   d,
#   theta[patch, year],
#   alpha,
#   init_dep[patch],
#   sigma_obs,
#   sigma_r_raw,
#   log_r0,
#   p_length_50_sel,
#   Topt,
#   width,
#   beta_obs_int,
#   beta_obs
# ) 
diagnostic_frame <- test |> 
  left_join(diagnostics) |> 
  mutate(divergent = divergent__) |> 
  mutate(rando = rnorm(length(divergent)))

diagnostics |> 
  group_by(.chain) |> 
  summarise(d = mean(divergent__ == 1))

diagnostic_frame |> 
  ggplot(aes(sigma_obs, log_r0, color = divergent__)) +
  geom_point()

diagnostic_frame |> 
  ggplot(aes(log_r0, divergent__)) +
  geom_point()

# divergent_drivers <- ranger(divergent ~ ., data = diagnostic_frame |> select(-contains("__"),-starts_with(".")), importance = "permutation", num.trees = 5000, classification = TRUE)

# vip(divergent_drivers,
#     num_features = 20)
# library(rattle)
# library(rpart.plot)
# library(RColorBrewer)

# plot mytree
# divergent_drivers <- rpart(
#   divergent ~ .,
#   data = diagnostic_frame |> select(-contains("__")),
#   method = "class",
#   minsplit = 2,
#   minbucket = 10000
# )

# fancyRpartPlot(divergent_drivers, caption = NULL)

test <- tidybayes::spread_draws(diagnostic_fit, density_hat[patch, year], theta[patch, year], sigma_obs)

habitat <- tidybayes::spread_draws(diagnostic_fit, habitat[patch])

habitat |>
  ggplot(aes(habitat)) +
  geom_histogram() +
  facet_wrap(~patch, scales = "free")

load(here("processed-data","stan_data_prep.Rdata"))



# composition data

length_comps <- reshape2::melt(len) |>
  rename(
    patch = Var1,
    length_bin = Var2,
    year = Var3,
    count = value
  ) |>
  as_tibble() |>
  group_by(year, patch) |>
  mutate(pcount = count / sum(count))  |>
  mutate(pcount = replace_na(pcount, 0)) |> 
  ungroup()

length_comps |> 
  filter(year == max(year)) |> 
  ggplot(aes(length_bin, pcount)) +
  geom_point() +
  facet_wrap(~patch)


tmp <- tidybayes::spread_draws(diagnostic_fit, n_at_length_hat[year, patch, length_bin], ndraws = 200)

  
estimated_length_comps <- tmp |>
  group_by(patch, year, .chain, .draw) |>
  mutate(pcount_hat = n_at_length_hat / sum(n_at_length_hat))

estimated_length_comps |> 
  filter(year == max(year)) |> 
  ggplot(aes(length_bin, pcount_hat)) + 
  stat_lineribbon(alpha = 0.25) +
  geom_point(data = length_comps |> filter(year == max(year)), aes(length_bin, pcount)) +
  facet_wrap(~patch , scales = "free_y") + 
  scale_fill_brewer()


# abundance trends --------------------------------------------------------

abund_p_y <-  dens %>%
  as.data.frame() |>
  mutate(patch = 1:np) |>
  pivot_longer(
    -patch,
    names_to = "year",
    values_to = "abundance",
    names_prefix = "V",
    names_transform = list(year = as.integer)
  )


sbt_p_y <-  sbt %>%
  as.data.frame() |>
  mutate(patch = 1:np) |>
  pivot_longer(
    -patch,
    names_to = "year",
    values_to = "sbt",
    names_prefix = "V",
    names_transform = list(year = as.integer)
  )

abund_p_y <- abund_p_y |> 
  left_join(sbt_p_y, by = c("year", "patch"))

abund_p_y |> 
  ggplot(aes(sbt, abundance)) + 
  geom_point() + 
  facet_wrap(~patch)


abund_p_y |> 
  ggplot(aes(year, abundance)) + 
  geom_point() + 
  facet_wrap(~patch)

abund_p_y |> 
  group_by(patch) |> 
  summarise(any_not_zero = any(abundance > 0))

abund_p_y |> 
  ggplot(aes(year, abundance == 0)) + 
  geom_point() + 
  facet_wrap(~patch)

abund_p_y |> 
  group_by(year) |> 
  summarise(abundance = sum(abundance)) |> 
  ggplot(aes(year, abundance)) + 
  geom_point() +
  scale_y_continuous(limits = c(0, NA))


together <- test |> 
  mutate(predicted_abundance = density_hat / theta) |> 
  left_join(abund_p_y, by = c("patch", "year"))


together |> 
  ggplot(aes((abundance), (predicted_abundance))) + 
  geom_point(alpha = 0.25) + 
  geom_abline(slope = 1, intercept = 0, color = "red") + 
  geom_smooth(method = "lm") + 
  scale_x_continuous("Observed Number Density") + 
  scale_y_continuous(name = "Probability of occurance * Predicted Number Density") + 
  theme_minimal()


abund_p_y_hat <- tidybayes::spread_draws(diagnostic_fit, density_hat[patch,year], ndraws  =100)

abundance_v_time <- abund_p_y_hat %>%
  ggplot(aes(year, density_hat)) +
  stat_lineribbon() +
  geom_point(data = abund_p_y, aes(year, abundance), color = "red") +
  facet_wrap(~patch) +
  labs(x="Year",y="Abundance", fill="Probability") +
  scale_fill_brewer()

abundance_v_time

