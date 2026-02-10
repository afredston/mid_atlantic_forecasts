#################
# This script assembles all of the parameters needed to run the DRM for summer flounder 
# It is hard-coded by nature and would need to be substantially revised for new species 

set.seed(42)
library(tidyverse)
library(here)
library(brms) # used for interpolating a few temperature values
funs <- list.files("functions")
sapply(funs, function(x) source(file.path("functions",x)))

dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
hauldat <- dat |> 
  select(haulid, btemp, year, date, lat, lon) |> 
  distinct()
dat <- dat %>% filter(length >17)

dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
hauldat_test <- dat_test |> 
  select(haulid, btemp, year, date, lat, lon) |> 
  distinct()

dat_f_age_prep <- read_csv(here("processed-data","summer_flounder_F_by_age.csv")) %>%
  rename_with(str_to_lower)
wt_at_age_raw <- read_csv(here("processed-data","summer_flounder_wt_at_age.csv"))
#############
# make model decisions that involve data prep
#############

# set fixed parameters from stock assessment (NOAA SAW 66)
loo = 83.6
k = 0.14
m = 0.25
age_at_maturity = 3 # https://www.fisheries.noaa.gov/species/summer-flounder
t0=-.2
cv= 0.2 # guess
min_age = 0
max_age = 15 
length_50_sel_guess= 20 
age_sel= 0
sel_100 = 3 # this is actually age 2, but we start counting ages at 0 (recruits), hence the 3 passed to Stan here 
h = 0.8

# the f-at-age data starts in 1982; fill in the previous years with the earliest year of data
f_early <- expand_grid(year=seq(1972, 1981, 1), age=unique(dat_f_age_prep$age)) %>% 
  left_join(dat_f_age_prep %>% filter(year==1982) %>% select(age, f)) 
dat_f_age_prep <- bind_rows(dat_f_age_prep, f_early)

##########
# prep data for fitting
##########

# reshape fish data 
patches <- sort(unique(dat %>% 
                         mutate(lat_floor = floor(lat), .keep="none") %>% 
                         pull(lat_floor)))

np = length(patches) 

dat_train_lengths <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% patches)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_lengths <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  group_by(length, year, lat_floor) %>% 
  summarise(sum_num_at_length = sum(number_at_length)) %>% 
  filter(lat_floor %in% patches)%>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_train_dens <- dat %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% patches) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size (often zeros) 
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo, including zeros 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_dens <- dat_test %>% 
  mutate(lat_floor = floor(lat)) %>% 
  filter(lat_floor %in% patches) %>% 
  group_by(haulid) %>% 
  mutate(dens = sum(number_at_length)) %>% # get total no. fish in each haul, of any size
  group_by(year, lat_floor) %>% 
  summarise(mean_dens = mean(dens)) %>%  # get mean density (all sizes) / haul for the patch*year combo 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

# get time dimension
years <- sort(unique(dat_train_lengths$year)) 
years_proj <- sort(unique(dat_test_lengths$year))
ny <- length(years)
ny_proj <- length(years_proj)

# get temperature data

# how many data points per patch and year? 
nhauls_dat_sbt <- hauldat |> 
  filter(!is.na(btemp)) |> 
  mutate(lat_floor = floor(lat)) |> 
  group_by(year, lat_floor) |> 
  summarise(n=n())

nhauls_dat_test_sbt <- hauldat_test |> 
  filter(!is.na(btemp)) |> 
  mutate(lat_floor = floor(lat)) |> 
  group_by(year, lat_floor) |> 
  summarise(n=n())

quantile(nhauls_dat_sbt$n, probs=0.05)
quantile(nhauls_dat_test_sbt$n, probs = 0.05)
# drop the lowest 5% of data, i.e., when there are <= 3 records per patch 

sbt_patches_dat_ok <- nhauls_dat_sbt |> 
  filter(n>3) |> 
  mutate(key = paste0(year,"-",lat_floor))
sbt_patches_dat_test_ok <- nhauls_dat_test_sbt |> 
  filter(n>3)|> 
  mutate(key = paste0(year,"-",lat_floor))

# get temperature data if there is enough data ... 
dat_train_sbt <- hauldat |> 
  mutate(lat_floor = floor(lat),
         key = paste0(year,"-",lat_floor)) %>% 
  filter(key %in% sbt_patches_dat_ok$key) |> 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

dat_test_sbt <- hauldat_test |> 
  mutate(lat_floor = floor(lat),
         key = paste0(year,"-",lat_floor)) %>% 
  filter(key %in% sbt_patches_dat_test_ok$key) |> 
  group_by(lat_floor, year) %>% 
  summarise(sbt = mean(btemp, na.rm=TRUE)) %>% 
  ungroup() %>% 
  mutate(patch = as.integer(as.factor(lat_floor)))

# check that dataframes are complete--they are not because some temperature data are missing 
nrow(dat_train_sbt)==(np*ny) # false because of missing data
nrow(dat_test_sbt) == (np*ny_proj) # false because of missing data

# use one model to predict into both 
sbt_brm <- brm(
  btemp ~ lat + (1 | year), # note that this leverages the raw latitude and temperature data, not aggregated to patches 
  data   = bind_rows(hauldat, hauldat_test),
  family = gaussian(),
  cores  = 4
) 

keys_train <- expand_grid(year = years, lat_floor = patches) |> 
  mutate(key = paste0(year,"-",lat_floor)) 
keys_test <- expand_grid(year = years_proj, lat_floor = patches) |> 
  mutate(key = paste0(year,"-",lat_floor)) 

# make dataframes of missing data to fill in; this will also identify year*patch combos with zero data, as well as those that don't meet the threshold set above for number of data points
tmp_fill_train <- keys_train |> 
  filter(!key %in% sbt_patches_dat_ok$key) |> 
  mutate(lat = lat_floor + 0.5) # predict into center of patch 
# note that these are mostly patches 1 and 10 missing 

tmp_fill_test <- keys_test |> 
  filter(!key %in% sbt_patches_dat_test_ok$key) |> 
  mutate(lat = lat_floor + 0.5)  # predict into center of patch 
# patch 10 often missing but also a data gap in 2007/2008

# predict missing data
pred_fill_train <- fitted(
  sbt_brm, 
  newdata = tmp_fill_train,
  re_formula = NULL,
  summary = TRUE
)

pred_fill_test <- fitted(
  sbt_brm, 
  newdata = tmp_fill_test,
  re_formula = NULL,
  summary = TRUE
)

# combine with lat/lon/year info 
sbt_fill_train <- cbind(tmp_fill_train, pred_fill_train) |> 
  select(year, lat_floor, Estimate) |> 
  rename(sbt = Estimate) |> 
  mutate(patch = lat_floor - min(patches) + 1)
sbt_fill_test <- cbind(tmp_fill_test, pred_fill_test) |> 
  select(year, lat_floor, Estimate) |> 
  rename(sbt = Estimate)|> 
  mutate(patch = lat_floor - min(patches) + 1)

dat_test_sbt <- bind_rows(dat_test_sbt, sbt_fill_test)
nrow(dat_test_sbt) == (np*ny_proj) # should be true now 

dat_train_sbt <- bind_rows(dat_train_sbt, sbt_fill_train)
nrow(dat_train_sbt) == (np*ny) # should be true now 

# make length to age conversions
length_at_age_key <-
  generate_length_at_age_key(
    min_age = min_age,
    max_age = max_age,
    cv = cv,
    linf = loo,
    k = k,
    t0 = t0,
    time_step = 1,
    linf_buffer = 1.5
  )

l_at_a_mat <- length_at_age_key %>% 
  select(age, length_bin, p_bin) %>% 
  pivot_wider(names_from = length_bin, values_from = p_bin) %>% 
  ungroup() %>% 
  select(-age) %>% 
  as.matrix()

# tidy f data
# the source data only has f estimates up to age 7
# fill in other ages with f
older_ages <- expand_grid(age=seq(max(dat_f_age_prep$age)+1, max_age, 1), year= unique(dat_f_age_prep$year)) %>% 
  left_join(dat_f_age_prep %>% filter(age==max(age)) %>% select(year, f))
dat_f_age_prep <- dat_f_age_prep %>% 
  bind_rows(older_ages)

dat_f_age <- dat_f_age_prep %>% 
  filter(year %in% years) 

dat_f_age_proj <- dat_f_age_prep %>% 
  filter(year %in% years_proj) %>%
  bind_rows(dat_f_age %>% filter(year==max(year))) # need final year of training data to initialize projection

lbins <- unique(length_at_age_key$length_bin)
n_lbins <- length(lbins) 

n_ages <- nrow(l_at_a_mat)

# now that we have n_ages, calculate weight at age
wt_at_age_prep <- wt_at_age_raw %>% 
  filter(!Age %in% seq(7, 10, 1)) %>% 
  mutate(Age = gsub("over7",7,Age),
         Age = as.numeric(Age),
         Age = Age + 1) %>% # start at 1 not 0
  group_by(Age) %>% 
  summarise(wt = mean(Wt)) %>% # average over all years 
  ungroup() %>% 
  arrange(Age)
# fill in plus ages
wt_at_age_add <- data.frame(Age = seq(max(wt_at_age_prep+1),n_ages,1), wt = slice_tail(wt_at_age_prep)$wt)
wt_at_age <- rbind(wt_at_age_prep, wt_at_age_add)$wt

# now that years are defined above, convert them into indices in the datasets
# be sure all these dataframes have exactly the same year range! 

dat_train_dens$year = as.integer(as.factor(dat_train_dens$year))
dat_test_dens$year = as.integer(as.factor(dat_test_dens$year))
dat_train_lengths$year = as.integer(as.factor(dat_train_lengths$year))
dat_test_lengths$year = as.integer(as.factor(dat_test_lengths$year))
dat_test_sbt$year= as.integer(as.factor(dat_test_sbt$year))
dat_train_sbt$year= as.integer(as.factor(dat_train_sbt$year))
dat_f_age$year = as.integer(as.factor(dat_f_age$year))
dat_f_age_proj$year = as.integer(as.factor(dat_f_age_proj$year))

# make matrices/arrays from dfs -- slow! 
len <- array(0, dim = c(np, n_lbins, ny)) 
for(p in 1:np){
  for(l in 1:n_lbins){
    for(y in 1:ny){
      tmp <- dat_train_lengths %>% filter(patch==p, round(length)==lbins[l], year==y) 
      if (nrow(tmp) > 0){
        len[p,l,y] <- tmp$sum_num_at_length
      }
    }
  }
}

# plot(len[4,,20])

dens <- array(NA, dim=c(np, ny))

  for(p in 1:np){
    for(y in 1:ny){
      tmp2 <- dat_train_dens %>% filter(patch==p, year==y) 
      dens[p,y] <- tmp2$mean_dens #  * (1/0.0384) * patcharea
      # previously converting fish/tow to fish, mean counts (in fish/tow) * tows/km2 * km2 ==> fish 
    # now just modeling densities in the net 
      }
  }

sbt <- array(NA, dim=c(np,ny))
for(p in 1:np){
  for(y in 1:ny){
    tmp3 <- dat_train_sbt %>% filter(patch==p, year==y) 
    sbt[p,y] <- tmp3$sbt
  }
}

sbt_proj <- array(NA, dim=c(np,ny_proj))
for(p in 1:np){
  for(y in 1:ny_proj){
    tmp6 <- dat_test_sbt %>% filter(patch==p, year==y) 
    sbt_proj[p,y] <- tmp6$sbt
  }
}


f <- array(NA, dim=c(n_ages,ny))
for(a in min_age:max_age){
  for(y in 1:ny){
    tmp4 <- dat_f_age %>% filter(age==a, year==y) 
    f[a+1,y] <- tmp4$f # add 1 because matrix indexing starts at 1 not 0
  } 
}

f_proj <- array(NA, dim=c(n_ages,(ny_proj+1)))
for(a in min_age:max_age){
  for(y in 1:(ny_proj+1)){
    tmp5 <- dat_f_age_proj %>% filter(age==a, year==y) 
    f_proj[a+1,y] <- tmp5$f # add 1 because matrix indexing starts at 1 not 0
    
  }
}

a <- seq(min_age, max_age)

check <- a %*% l_at_a_mat

bin_mids=lbins+0.5 

save(
  dat_train_dens,
  dat_test_dens,
  np,
  n_ages,
  ny,
  ny_proj,
  n_lbins,
  len,
  dens,
  sbt,
  sbt_proj,
  m,
  f,
  f_proj,
  k,
  loo,
  h,
  t0,
  cv,
  length_50_sel_guess,
  n_lbins, 
  age_sel,
  bin_mids,
  sel_100=sel_100,
  age_at_maturity,
  l_at_a_mat,
  wt_at_age,
  patches, 
  years,
  years_proj,
  file=here("processed-data","stan_data_prep.Rdata")
)
