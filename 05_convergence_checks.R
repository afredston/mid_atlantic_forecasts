set.seed(42)
library(tidyverse)
#library(tidybayes)
library(here)
library(rstan)
#library(rstanarm)
library(cmdstanr)

ctrl_file <- read_csv("ctrl_file_used.csv") 

run_in_parallel <- TRUE

num_iters_per_chain <- 3000 # no longer part of the fitted stan object so copying over here from run_drms (iters - warmups)

if(run_in_parallel == TRUE){
  library(parallel)
  num_cores <- 40
}

# make diagnostic files 
if(run_in_parallel == TRUE){
write_out_diagnostics <- function(k) {
  i <- ctrl_file$id[k]
  results_path <- here("results", i)
  stanmod <- read_rds(file.path(results_path, "stan_model_fit.rds"))
  diagnostic_fit <- stanmod$diagnostic_summary()
  write_rds(diagnostic_fit, file = file.path(results_path, "diagnostic_fit.rds"))
  rm(stanmod)
}

cl <- makeCluster(num_cores) 
clusterExport(cl, c("ctrl_file", "write_out_diagnostics")) 

clusterEvalQ(cl, {
  library(readr)
  library(here)
  }) 

parLapply(cl, 1:nrow(ctrl_file), write_out_diagnostics) 

stopCluster(cl)

} else {
for(k in 1:nrow(ctrl_file)){
  i = ctrl_file$id[k]  
  results_path <- file.path("results",i)
  diagnostic_fit <- read_rds(paste0(results_path, "/stan_model_fit.rds"))
  diagnostic_fit <- diagnostic_fit$diagnostic_summary()
  write_rds(diagnostic_fit, file = paste0(results_path, "/diagnostic_fit.rds"))
}
} 

# summarize diagnostics and calculate convergence checks 

datalist = list()
for(k in 1:nrow(ctrl_file)){
  
  i = ctrl_file$id[k]  
  
  results_path <- here("results",i)
  
  diagnostic_ls <- readRDS(here(results_path,"diagnostic_fit.rds"))
  
  dat <- data.frame(id = i,
    successful_chains = length(diagnostic_ls$num_divergent),
                    mean_divergences = sum(diagnostic_ls$num_divergent) / (num_iters_per_chain*length(diagnostic_ls$num_divergent)),
                    low_ebfmi = ifelse(mean(diagnostic_ls$ebfmi)<0.3, TRUE, FALSE))
  
  datalist[[k]] <- dat
}

summarydat = do.call(rbind, datalist)

write_csv(summarydat, file=here("results","convergence_checks.csv"))
