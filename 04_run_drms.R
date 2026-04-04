# for compilation issues with cmdstanr try cmdstanr::rebuild_cmdstan()
set.seed(42)
library(tidyverse)
library(ggridges)
library(tidybayes)
library(here)
library(rstan)
library(Matrix)
library(cmdstanr)
library(data.table) 
library(parallel)
set.seed(424242)
rstan_options(javascript = FALSE, auto_write = TRUE)

run_in_parallel <- TRUE
if(run_in_parallel == TRUE){
  n_cores <- 160  # set cores 
}

# load fit_drm function 
# fit_drm() fits the model and writes out the model object and a plot to the results directory
funs <- list.files("functions")
sapply(funs, function(x)
  source(file.path("functions", x)))

# which range edges should be calculated?
#  
quantiles_calc <- c(0.05, 0.5, 0.95)

ctrl_file <- read_csv("control_file.csv") %>% 
  filter(process_error_toggle == 1, 
         eval_l_comps == 0,
         known_f == 1,
         spawner_recruit_relationship == 0)

ctrl_file$name <- c("DRM null", "DRM T-movement", "DRM T-recruit", "DRM T-mortality" )

write_csv(ctrl_file, file=here("ctrl_file_used.csv"))
# ctrl_file <- read_csv("control_file.csv") %>%
#   filter(
#     eval_l_comps == 0,
#     spawner_recruit_relationship == 1,
#     process_error_toggle == 1,
#     known_f == 0,
#     T_dep_mortality == 0
#   ) |>
#   ungroup() |>
#   slice(1)

fit_drms <- TRUE
use_poisson_link <- 0
# if (use_poisson_link){
#   run_name <- "yes-pois"
# } else {
#   run_name <- "no-pois"
# }
write_summary <- TRUE
iters <- 5000
warmups <- 2000
chains <- 4
cores <- 4

load(here("processed-data","stan_data_prep.Rdata"))

if(run_in_parallel == TRUE){
  run_drm_fn <- function(k) {
    i <- ctrl_file$id[k]
    results_path <- file.path("results", i)
    
    if (fit_drms == TRUE) {
      drm_fits <- ctrl_file %>%
        filter(id == i)
      
      drm_fits$fits <- list(tryCatch(fit_drm(
        use_poisson_link = use_poisson_link,
        create_dir = TRUE,
        run_name = drm_fits$id,
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
        quantiles_calc = quantiles_calc
      ))) 
    }
  }
  
  mclapply(1:nrow(ctrl_file), run_drm_fn, mc.cores = n_cores) # does not return anything just writes out models 
} else {
  for(k in 1:nrow(ctrl_file)){
    i = ctrl_file$id[k]  
    
    results_path <- file.path("results",i)
    
    # turn off if you just want to load already-fitted models and analyze them
    
    if (fit_drms==TRUE){
      drm_fits <-  ctrl_file %>%
        filter(id == i)
      
      drm_fits$fits <- list(tryCatch(fit_drm(
        use_poisson_link = use_poisson_link,
        create_dir = TRUE,
        run_name = drm_fits$id,
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
}
