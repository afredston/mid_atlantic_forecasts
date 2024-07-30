
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

ctrl_file <- read_csv("control_file.csv")
ctrl_file <- read_csv("control_file.csv") %>%
  filter(
    eval_l_comps == 0,
    spawner_recruit_relationship == 0,
    process_error_toggle == 0,
    known_f == 1,
    T_dep_mortality == 1
  ) |>
  ungroup() |>
  slice(1)

fit_drms <- TRUE
use_poisson_link <- TRUE
if (use_poisson_link){
  run_name <- "yes-pois"
} else {
  run_name <- "no-pois"
}
make_plots <- TRUE
write_summary <- TRUE
iters <- 10000
warmups <- 5000
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
      adapt_delta = 0.99, 
      run_forecast = 1,
      quantiles_calc = quantiles_calc, 
    )
    ) 
    )# as currently written this just adds a column to drm_fits that says "all done". the column `fits` used to contain the model object itself 
    
    
  } # close fit_drms 
}

diagnostic_fit <- read_rds(here("results",run_name, "stan_model_fit.rds"))

diagnostic_fit$diagnostic_summary()

test <- tidybayes::spread_draws(diagnostic_fit, density_hat[patch,year],theta[patch,year], ndraws  =100)

habitat <- tidybayes::spread_draws(diagnostic_fit, habitat[patch], ndraws  =100)

habitat |> 
  ggplot(aes(habitat)) + 
  geom_histogram() +
  facet_wrap(~patch, scales = "free")

load(here("processed-data","stan_data_prep.Rdata"))


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


  together <- test |> 
    mutate(predicted_abundance = density_hat / theta) |> 
    left_join(abund_p_y, by = c("patch", "year"))
  
  
  together |> 
    ggplot(aes(abundance, predicted_abundance)) + 
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

