
ctrl_file <- read_csv(file=here("ctrl_file_used.csv"))
run_in_parallel <- TRUE

if(run_in_parallel == TRUE){
  library(parallel)
  num_cores <- 40
}

# thin and tidy posteriors 
if(run_in_parallel == TRUE){
  thin_posteriors <- function(k) {
    i <- ctrl_file$id[k]
    results_path <- here("results", i)
    tmp_model <- read_rds(file.path(results_path, "stan_model_fit.rds"))
    
    # first get the static posteriors 
    d <- gather_draws(tmp_model, d) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(.value)) %>%  # average across chains
      select(value) %>% 
      mutate(param = "d")
    
    Topt <- gather_draws(tmp_model, Topt) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(.value)) %>%  # average across chains
      select(value) %>% 
      mutate(param = "Topt")
    
    width <- gather_draws(tmp_model, width) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(.value)) %>%  # average across chains
      select(value) %>% 
      mutate(param = "width")
    
    sigma_obs <- gather_draws(tmp_model, sigma_obs) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(.value)) %>%  
      select(value) %>% 
      mutate(param = "sigma_obs")
    
    mean_recruits <- gather_draws(tmp_model, log_mean_recruits) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(.value)) %>%  
      select(value) %>% 
      mutate(param = "mean_recruits", value = exp(value)) # exponentiate so not in log space 
    
    tmp <- rbind(d, Topt, width, sigma_obs, mean_recruits) 
    write_rds(tmp, file = file.path(results_path, "fixed_params_averaged.rds"))
    
    # now thin the latent states (pop dy over space and time) 
    density_hat_thin <- read_rds(here(results_path, "density_hat.rds"))  |> 
      filter(.draw %% 50 == 0) # thin to 1 sample every 50 iterations on each chain 
    write_rds(density_hat_thin, file = file.path(results_path, "density_hat_thin.rds"))
    
    centroid_proj_thin <- read_rds(here(results_path, "centroid_proj.rds"))  |> 
      filter(.draw %% 50 == 0) # thin to 1 sample every 50 iterations on each chain 
    write_rds(centroid_proj_thin, file = file.path(results_path, "centroid_proj_thin.rds"))
    
    density_obs_proj_thin <- read_rds(here(results_path, "density_obs_proj.rds"))  |> 
      filter(.draw %% 50 == 0) # thin to 1 sample every 50 iterations on each chain 
    write_rds(density_obs_proj_thin, file = file.path(results_path, "density_obs_proj_thin.rds"))
    
    range_quantiles_proj_thin <- read_rds(here(results_path, "range_quantiles_proj.rds"))  |> 
      filter(.draw %% 50 == 0) # thin to 1 sample every 50 iterations on each chain 
    write_rds(range_quantiles_proj_thin, file = file.path(results_path, "range_quantiles_proj_thin.rds"))
    
    range_quantiles_thin <- read_rds(here(results_path, "range_quantiles.rds"))  |> 
      filter(.draw %% 50 == 0) # thin to 1 sample every 50 iterations on each chain 
    write_rds(range_quantiles_thin, file = file.path(results_path, "range_quantiles_thin.rds"))
    
    sigma_r <- tidybayes::spread_draws(tmp_model, sigma_r[year]) %>% 
      group_by(.iteration) %>% 
      summarise(value = mean(sigma_r)) %>%  
      select(value) %>% 
      mutate(param = "sigma_r")
    
    write_rds(sigma_r, file = file.path(results_path, "sigma_r.rds"))
    
    rm(list=ls())
  }
  
  cl <- makeCluster(num_cores) 
  clusterExport(cl, c("ctrl_file", "thin_posteriors")) 
  
  clusterEvalQ(cl, {
    library(tidyverse)
    library(here)
    library(tidybayes)
  }) 
  
  parLapply(cl, 1:nrow(ctrl_file), thin_posteriors) 
  
  stopCluster(cl)
  
} 

