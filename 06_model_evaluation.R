###############
# SETUP
###############

# load packages 
set.seed(42)
library(tidyverse)
library(magrittr)
library(here)
library(mgcv)
library(tidybayes)
library(ggrepel) # move plots to another script 
#library(parallel)
library(doParallel)
# read in data 
convergence_checks <- read_csv(file=here("results","convergence_checks.csv"))
ctrl_file <- read_csv(file=here("ctrl_file_used.csv"))
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
load(here("processed-data","stan_data_prep.Rdata"))
quantiles_calc <- c(0.05, 0.5, 0.95)
run_in_parallel <- TRUE
if(run_in_parallel == TRUE) {  num_cores <- 40 }

# haul_km2_per_tow <- 0.0384 # conversion factor from hauls to km2; this is in units of km2/tow. 
# see: https://github.com/afredston/marine_heatwaves_trawl/blob/main/prep_trawl_data.R#L59 


# make user decisions 
divergence_cutoff <- 0.05
chains_cutoff <- 3 
generate_exploratory_plots <- FALSE

# identify models that pass convergence checks 
summarydat <- convergence_checks %>% 
  left_join(ctrl_file) %>% 
  filter(mean_divergences <= divergence_cutoff, 
         successful_chains >= chains_cutoff)

# which models did not pass convergence checks?
ctrl_file %>% 
  filter(!id %in% summarydat$id) %>% 
  pull(description, id)

# write a function to calculate patch edge positions 
calculate_range_edge <- function(patches, weights, q){
  if(length(patches) == length(weights)){
    cutoff <- q * sum(weights)
    csum <- cumsum(weights)
    p <- min(which(csum > cutoff)) 
    if(p == 1) { # special case when range edge is in the first patch 
      out <- patches[1] + cutoff / weights[1]
    } else {
      dec <- (cutoff - csum[p-1]) / weights[p] # calculate decimal position into the edge patch -- what proportion of the patch would be occupied, assuming a constant density 
      out <- patches[p] + dec 
      return(out) }
  } else {return(print("Length of patches and weights must be equal!"))}
}

# calculate summary statistics on real data 
dat_train_patch <- dat_train_dens %>% 
  group_by(year) %>% 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.05),
    centroid = weighted.mean(lat_floor, w=mean_dens),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.95),
    dens = mean(mean_dens)) %>% 
  arrange(year) %>% 
  mutate(year = year + min(years) - 1, 
         dens_lr = log(dens / lag(dens))) %>% 
  pivot_longer(cols=c(warm_edge:dens_lr), names_to="feature", values_to="value") 
write_csv(dat_train_patch, file=here("processed-data","dat_train_patch.csv"))


dat_test_patch <- dat_test_dens %>% 
  group_by(year) %>% 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.05),
    centroid = weighted.mean(lat_floor, w=mean_dens),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.95),
    dens = mean(mean_dens)) %>% 
  arrange(year) %>% 
  mutate(year = year + min(years_proj) - 1, 
         dens_lr = log(dens / lag(dens))) %>% 
  pivot_longer(cols=c(warm_edge:dens_lr), names_to="feature", values_to="value") 

# note that one dens_lr value is missing; let's calculate it manually here 
dat_test_patch[dat_test_patch$year==min(years_proj) & dat_test_patch$feature=="dens_lr",]$value <- log(dat_test_patch[dat_test_patch$year==min(years_proj) & dat_test_patch$feature=="dens",]$value / dat_train_patch[dat_train_patch$year==max(years) & dat_train_patch$feature=="dens",]$value)

write_csv(dat_test_patch, file=here("processed-data","dat_test_patch.csv"))

# create other df for the entire time-series 
time_series_dat <- dat_test_dens %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = lat_floor, Density=mean_dens, .keep="none") %>%
  bind_rows(dat_train_dens |> mutate(Year = (year + min(years) - 1), Latitude = lat_floor, Density=mean_dens, .keep="none")) |> 
  group_by(Year) %>% 
  summarise(
    `Warm edge` = calculate_range_edge(patches=Latitude, weights=Density, q=0.05),
    Centroid = weighted.mean(Latitude, w=Density),
    `Cold edge` = calculate_range_edge(patches=Latitude, weights=Density, q=0.95)) %>% 
  arrange(Year) %>% 
  pivot_longer(cols=c(`Warm edge`:`Cold edge`), names_to="feature", values_to="value") 
write_csv(time_series_dat, file=here("processed-data","time_series_summary_stats.csv"))

###############
# FIT SDMS FOR COMPARISON
###############

# prep data for fitting GAM 

# fit GAM to all hauls (not patch level)
dat_train_gam <- dat %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length))

dat_test_gam <- dat_test %>% 
  group_by(haulid, year, lat, btemp) %>% 
  summarise(dens = sum(number_at_length))

# check data coverage of testing data
if(generate_exploratory_plots==TRUE){
  dat_test_gam %>% 
    mutate(lat_floor = floor(lat)) %>% 
    group_by(year, lat_floor) %>% 
    mutate(n_tot = n(), 
           test = is.na(btemp)) %>% 
    filter(test==FALSE) %>% 
    mutate(n_pres = n(), 
           prop_btemp = n_pres / n_tot) %>% 
    select(prop_btemp, year, lat_floor) %>% 
    distinct() %>% 
    ggplot() + 
    geom_tile(aes(x=year, y=lat_floor, fill=prop_btemp, color=prop_btemp)) }
# because so much environmental data is missing in 2008, drop it from predictions 

dat_test_gam %<>% 
  filter(!year==2008)

# fit GAMs 
# see https://github.com/pinskylab/project_velocity/blob/master/6_model_fitting_loop.R#L182

spdata <- dat_train_gam %>%
  mutate(pres = ifelse(dens>0, 1, 0), 
         logdens = log(dens)) %>% 
  filter(!is.na(btemp))
spdata_proj <- dat_test_gam %>%
  mutate(pres = ifelse(dens>0, 1, 0), 
         logdens = log(dens)) %>% 
  filter(!is.na(btemp))

mypresmodtt<-formula(pres ~ s(btemp))
myabunmodtt<-formula(logdens ~ s(btemp))

gammaPA <- log(nrow(spdata)) / 2
gammaAbun <- log(nrow(spdata[spdata$pres==1,])) / 2

mygam1tt<-gam(mypresmodtt, family="binomial",data=spdata, select=TRUE, gamma=gammaPA) 
mygam2tt<-gam(myabunmodtt, data=spdata[spdata$pres==1,], select=TRUE, gamma=gammaAbun) 

preds1tt <- predict(mygam1tt, newdata = spdata_proj, type="response") 
preds2tt <- exp(predict(mygam2tt, newdata = spdata_proj, type='response'))

predstt <- preds1tt*preds2tt 
predstt[predstt<0] = 0
predstt[is.na(predstt)] = 0 # is this correct?

spdata_proj$predstt <- predstt

gam_out <- spdata_proj %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(lat_floor, year) %>% 
  summarise(dens_pred = mean(exp(predstt)))# aggregate to patch scale for comparison to DRM 
write_csv(gam_out, file = here("processed-data","gam_density_time.csv"))

# calculate residuals by feature (centroid, edges) of forecast 

gam_time <- gam_out %>% 
  group_by(year) %>% # calculate summary stats 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.05),
    centroid = weighted.mean(lat_floor, w=dens_pred),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=dens_pred, q=0.95),
    dens = mean(dens_pred)) %>% 
  arrange(year) %>% 
  mutate(dens_lr = log(dens / lag(dens))) %>% # note that this is technically the LR of change from 2007 to 2009 because the GAM doesn't use 2008
  pivot_longer(cols=warm_edge:dens_lr, names_to="feature", values_to="value_tmp") 

gam_summary <- gam_time %>% 
  left_join(dat_test_patch)%>% # compare  to true data 
  mutate(resid = value_tmp - value, 
         resid_sq = resid^2, 
         .keep = "unused", # drop all the columns used in calculations 
         name = "GAM") 

###############
# MAKE PERSISTENCE FORECAST
###############

persistence <- dat_train_dens %>% 
  filter(year == max(year)) %>% 
  summarise(
    warm_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.05),
    centroid = weighted.mean(lat_floor, w=mean_dens),
    cold_edge = calculate_range_edge(patches=lat_floor, weights=mean_dens, q=0.95),
    dens = mean(mean_dens)) 

persistence_dat <- data.frame(year = years_proj) %>% 
  bind_rows(persistence) %>% 
  fill(warm_edge, centroid, cold_edge, dens, .direction="up") %>% 
  filter(!is.na(year))

persistence_summary <- persistence_dat %>% 
  pivot_longer(cols=warm_edge:dens, names_to="feature", values_to="value_tmp") %>% 
  left_join(dat_test_patch)%>% # compare  to true data 
  mutate(resid = value_tmp - value, 
         resid_sq = resid^2, 
         .keep = "unused", # drop all the columns used in calculations 
         name = "Persistence") 

points_for_plot <- dat_test_patch %>% 
  mutate(name = 'Observed')%>% rename(value_tmp = value) %>% 
  bind_rows(gam_time %>% mutate(name = 'GAM') )  %>% 
  bind_rows(persistence_dat %>% pivot_longer(cols=c('warm_edge','cold_edge','dens','centroid'), values_to='value_tmp', names_to='feature') %>% mutate(name='Persistence')) %>% 
  filter(feature %in% c('warm_edge','cold_edge','centroid')) %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature)) 
write_csv(points_for_plot, file=here("processed-data","points_for_plot.csv"))

###############
# SUMMARIZE DRM OUTPUTS
###############

# generate summary stats for the successful DRMs

drm_out <- NULL

# pull in drm forecasts and calculate residuals by year and by feature (centroid, edges) 
# slow! try to speed this up at some point

if(run_in_parallel == TRUE) {
  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  # Parallel loop
  drm_out <- foreach(i = 1:nrow(ctrl_file), .combine = bind_rows, .packages = c("here", "tidybayes", "dplyr", "readr")) %dopar% {
    tmpdat <- ctrl_file[i, ]
    
    results_path <- here(paste0("results/", tmpdat$id))
    
    tmp_model <- tryCatch(read_rds(file.path(results_path, "stan_model_fit.rds")), error = function(e) return(NULL))
    
    density_obs_proj <- tidybayes::spread_draws(tmp_model, density_obs_proj[patch,year])
    
    centroid_proj <- tidybayes::spread_draws(tmp_model, centroid_proj[year]) 
    
    range_quantiles_proj <- tidybayes::spread_draws(tmp_model, range_quantiles_proj[quantile, year]) %>%
      mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
    
    density_hat <- tidybayes::spread_draws(tmp_model, density_hat[patch,year])
    
    range_quantiles <- tidybayes::spread_draws(tmp_model, range_quantiles[quantile, year]) %>%
      mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
    
    # save those posteriors
    write_rds(density_obs_proj, file.path(results_path, "density_obs_proj.rds"))
    write_rds(density_hat, file.path(results_path, "density_hat.rds"))
    write_rds(range_quantiles_proj, file.path(results_path, "range_quantiles_proj.rds"))
    write_rds(range_quantiles, file.path(results_path, "range_quantiles.rds"))
    
    write_rds(centroid_proj, file.path(results_path, "centroid_proj.rds")) 
    
    centroid_tmp <- centroid_proj %>% 
      mutate(year = year + min(years_proj) - 1) %>% 
      left_join(dat_test_patch %>% filter(feature == "centroid")) %>% 
      mutate(resid = centroid_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name) %>% 
      select(-centroid_proj, -value)
    
    warm_edge_tmp <- range_quantiles_proj %>% 
      filter(range_quantiles_proj < Inf, quantile == 0.05) %>%  
      mutate(
        year = year + min(years_proj) - 1, 
        range_quantiles_proj = range_quantiles_proj + min(patches)) %>% 
      left_join(dat_test_patch %>% filter(feature == "warm_edge")) %>% 
      mutate(resid = range_quantiles_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name) %>% 
      select(-range_quantiles_proj, -value, -quantile)
    
    cold_edge_tmp <- range_quantiles_proj %>% 
      filter(range_quantiles_proj < Inf, quantile == 0.95) %>%  
      mutate(
        year = year + min(years_proj) - 1, 
        range_quantiles_proj = range_quantiles_proj + min(patches)) %>% 
      left_join(dat_test_patch %>% filter(feature == "cold_edge")) %>% 
      mutate(resid = range_quantiles_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name) %>% 
      select(-range_quantiles_proj, -value, -quantile)
    
    boop <- range_quantiles_proj %>% 
      filter(range_quantiles_proj < Inf, quantile == 0.95) %>%  
      mutate(
        year = year + min(years_proj) - 1, 
        range_quantiles_proj = range_quantiles_proj + min(patches)) %>% 
      left_join(dat_test_patch %>% filter(feature == "cold_edge")) %>% 
      mutate(resid = range_quantiles_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name) %>% 
      select(-value, -quantile) |> 
      group_by(year) |> 
      summarise(medlat = median(range_quantiles_proj))
    
    bind_rows(centroid_tmp, warm_edge_tmp, cold_edge_tmp)
  }
  
  # Stop the cluster
  stopCluster(cl)
}  else {
  for(i in 1:nrow(ctrl_file)){ # should be summarydat not ctrl_file once I update the code to do convergence checks
    
    tmpdat <- ctrl_file[i,]
    
    results_path <- here("results",tmpdat$id)
    
    # get the Stan model and extract posteriors that we want for plots 
    tmp_model <-  tryCatch(read_rds(file.path(results_path, "stan_model_fit.rds")))
    
    density_obs_proj <- tidybayes::spread_draws(tmp_model, density_obs_proj[patch,year])
    
    centroid_proj <- tidybayes::spread_draws(tmp_model, centroid_proj[year]) 
    range_quantiles_proj <- tidybayes::spread_draws(tmp_model, range_quantiles_proj[quantile, year]) %>%
      mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
    density_hat <- tidybayes::spread_draws(tmp_model, density_hat[patch,year])
    
    range_quantiles <- tidybayes::spread_draws(tmp_model, range_quantiles[quantile, year]) %>%
      mutate(quantile = as.factor(quantiles_calc[quantile]), .keep="unused") 
    # save those posteriors
    write_rds(density_obs_proj, file.path(results_path, "density_obs_proj.rds"))
    write_rds(density_hat, file.path(results_path, "density_hat.rds"))
    
    write_rds(range_quantiles_proj, file.path(results_path, "range_quantiles_proj.rds"))
    write_rds(range_quantiles, file.path(results_path, "range_quantiles.rds"))
    
    write_rds(centroid_proj, file.path(results_path, "centroid_proj.rds")) 
    
    centroid_tmp <- centroid_proj %>% 
      mutate(year = year + min(years_proj) - 1) %>% 
      left_join(dat_test_patch %>% filter(feature=="centroid")) %>% 
      mutate(resid = centroid_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name) %>% 
      select(-centroid_proj, -value)
    
    warm_edge_tmp <- range_quantiles_proj %>% 
      filter(range_quantiles_proj < Inf,
             quantile == 0.05) %>%  
      mutate(
        year = year + min(years_proj) - 1, 
        range_quantiles_proj = range_quantiles_proj + min(patches)) %>% # don't need to subtract 1. we want patch 0.09 to be lat 35.09
      left_join(dat_test_patch %>% filter(feature=="warm_edge")) %>% 
      mutate(resid = range_quantiles_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name)%>% 
      select(-range_quantiles_proj, -value, -quantile)
    
    cold_edge_tmp <- range_quantiles_proj %>% 
      filter(range_quantiles_proj < Inf,
             quantile == 0.95) %>%  
      mutate(
        year = year + min(years_proj) - 1, 
        range_quantiles_proj = range_quantiles_proj + min(patches)) %>% # don't need to subtract 1. we want patch 0.09 to be lat 35.09
      left_join(dat_test_patch %>% filter(feature=="cold_edge")) %>% 
      mutate(resid = range_quantiles_proj - value,
             resid_sq = resid^2, 
             name = tmpdat$name)%>% 
      select(-range_quantiles_proj, -value, -quantile)
    
    drm_out <- bind_rows(drm_out, centroid_tmp, warm_edge_tmp, cold_edge_tmp)
  }
}
write_csv(drm_out, file=here("processed-data","posteriors_for_model_evaluation.csv"))

# drm_out <- drm_out %>% 
#   filter(id %in% summarydat$id)

drm_summary <- drm_out %>% 
  filter(!is.na(resid)) %>% 
  group_by(year, feature, name) %>% 
  summarise(resid = mean(resid)) %>% 
  mutate(resid_sq = resid^2)

dat_forecasts <- drm_summary %>% 
  left_join(ctrl_file %>% select(name)) %>% 
  bind_rows(gam_summary, persistence_summary)%>% 
  filter(feature %in% c('centroid','cold_edge','warm_edge'))  # remove density
write_csv(dat_forecasts, file=here("processed-data","model_performance_by_year.csv"))

# pool across years to calculate bias and RMSE
dat_forecasts_summ <- dat_forecasts %>% 
  group_by(feature, name) %>% 
  summarise(RMSE = sqrt(mean(resid_sq)), 
            Bias = mean(resid)) %>% 
  pivot_longer(cols=c(RMSE, Bias), values_to="value", names_to="metric") 

write_csv(dat_forecasts_summ, file = here("processed-data","model_comparison_summary.csv"))

