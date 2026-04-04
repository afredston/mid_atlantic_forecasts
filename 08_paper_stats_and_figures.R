# load packages and data 

library(tidyverse)
library(here)
# packages for map
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")

# packages for nice plots
library(ggplot2)
library(tidytext)
library(ggrepel)
library(ggh4x)
library(ggdist)
library(directlabels)
library(patchwork)
library(brms)
theme_set(theme_bw())

# data
dat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv"))
dat_test <- read_csv(here("processed-data","flounder_catch_at_length_fall_testing.csv"))
dat_catchonly <- read_csv(here("processed-data","flounder_catch_fall_training.csv"))
dat_test_catchonly <- read_csv(here("processed-data","flounder_catch_fall_testing.csv"))
convergence_checks <- read_csv(file=here("results","convergence_checks.csv"))
dat_forecasts_summ <- read_csv(file = here("processed-data","model_comparison_summary.csv"))
ctrl_file <- read_csv(file=here("ctrl_file_used.csv"))
dat_test_patch <- read_csv(file=here("processed-data","dat_test_patch.csv"))
gam_out <- read_csv(file = here("processed-data","gam_density_time.csv"))
points_for_plot <- read_csv(file=here("processed-data","points_for_plot.csv"))
time_series_dat <- read_csv(here("processed-data","time_series_summary_stats.csv"))
load(here("processed-data","stan_data_prep.Rdata"))

hauldat <- bind_rows(dat, dat_test) |> 
  select(date, lat, lon) |> 
  distinct()

fixed_param_dat <- NULL
for(i in ctrl_file$id){
  results_path <- here('results',i)
  out <- read_rds(file.path(here("results"),i,"fixed_params_averaged.rds"))
  out$name <- ctrl_file[ctrl_file$id==i,]$name
  fixed_param_dat <- rbind(fixed_param_dat, out)
}
drm_out <- read_csv(here("processed-data","posteriors_for_model_evaluation.csv")) 

run_in_parallel <- TRUE
if(run_in_parallel == TRUE){
  n_cores <- 100  # set cores 
}
############
# calculate all the statistics reported in-text (nothing gets saved in this code block, just printed to the console)
############

# did temperature change in the trawl surveys across all time? 
brm(btemp ~ year, 
    data = bind_rows(dat_catchonly, dat_test_catchonly) %>%   select(btemp, year, lat),
    family = gaussian(), 
    cores = 4, 
    chains = 4, 
    iter = 2000
)|>
  posterior::as_draws_df() |>
  (\(x) tibble::tibble(
    term = "year",
    median = median(x$b_year),
    hpdi_lower = coda::HPDinterval(coda::as.mcmc(x$b_year), prob = 0.95)[1],
    hpdi_upper = coda::HPDinterval(coda::as.mcmc(x$b_year), prob = 0.95)[2]
  ))()

# how many hauls? 
length(unique(dat$haulid)) + length(unique(dat_test$haulid)) #12203
nrow(dat_catchonly) + nrow(dat_test_catchonly) #12203 — should be the same as above

# how many positive encounters? 
nrow(dat_catchonly[dat_catchonly$abundance>0,]) + nrow(dat_test_catchonly[dat_test_catchonly$abundance>0,]) # 2308 

# what percentage of hauls were positive encounters? 
(nrow(dat_catchonly[dat_catchonly$abundance>0,]) + nrow(dat_test_catchonly[dat_test_catchonly$abundance>0,])) / (nrow(dat_catchonly) + nrow(dat_test_catchonly)) # 0.1891338

# how frequently were more than 10 individuals caught? 
nrow(bind_rows(dat_catchonly, dat_test_catchonly) %>% 
       filter(abundance >= 10)) / nrow(bind_rows(dat_catchonly, dat_test_catchonly)) # 0.034

# how many individuals caught? 
sum(dat_catchonly$abundance) # 7713
sum(dat_test_catchonly$abundance) # 6312

# how many individuals caught per year?
bind_rows(dat_catchonly, dat_test_catchonly) %>% 
  group_by(year) %>% 
  summarise(n=sum(abundance))

# how many models passed convergence checks?
divergence_cutoff <- 0.05
chains_cutoff <- 3 
nrow(convergence_checks %>% 
       filter(mean_divergences <= divergence_cutoff, 
              successful_chains >= chains_cutoff)) # all of them


# did summer flounder shift north in the testing interval? 
brm(value_tmp ~ year, 
    data = points_for_plot %>% 
      filter(name == 'Observed') %>% 
      filter(feature == 'Centroid'),
    family = gaussian(),
    cores = 4) |> 
  posterior::as_draws_df() |>
  (\(x) coda::HPDinterval(coda::as.mcmc(x$b_year), prob = 0.95))()

brm(value_tmp ~ year, 
    data = points_for_plot %>% 
      filter(name == 'Observed') %>% 
      filter(feature == 'Warm Edge'),
    family = gaussian(),
    cores = 4) |> 
  posterior::as_draws_df() |>
  (\(x) coda::HPDinterval(coda::as.mcmc(x$b_year), prob = 0.95))()

brm(value_tmp ~ year, 
    data = points_for_plot %>% 
      filter(name == 'Observed') %>% 
      filter(feature == 'Cold Edge'),
    family = gaussian(),
    cores = 4) |> 
  posterior::as_draws_df() |>
  (\(x) coda::HPDinterval(coda::as.mcmc(x$b_year), prob = 0.95))()


print(fixed_param_dat |> group_by(name, param) |> 
        summarise(
          median = median(value),
          hpdi_lower = coda::HPDinterval(coda::as.mcmc(value), prob = 0.95)[1],
          hpdi_upper = coda::HPDinterval(coda::as.mcmc(value), prob = 0.95)[2],
          .groups = "drop"
        ), n=28)

############
# make area map figure
############

states <- ne_states(country = c("United States of America","Canada"), returnclass = "sf")
select <- dplyr::select

load(here("processed-data","stan_data_prep.Rdata"))
recdat <- read_csv(here("processed-data","flounder_catch_at_length_fall_training.csv")) %>% 
  mutate(lat_floor = floor(lat)) %>% 
  select(lat_floor, lat, lon) %>% 
  group_by(lat_floor) %>% 
  summarise(minlon = min(lon),
            maxlon = max(lon))

recs <- data_frame(ymin = patches, 
                   ymax = patches + 0.99, 
                   xmin = recdat$minlon + 0.1, 
                   xmax = xmin + 2,
                   x = xmin + (xmax-xmin)/2, 
                   y = ymin + 0.5
) # be sure to note that rectangle width is approximate

recs2 <- recs %>% 
  select(y, xmax, xmin)  %>%
  pivot_longer(cols = c(xmax, xmin), names_to="identity", values_to="x") %>% 
  group_by(y) %>% 
  mutate(meanx = mean(x),
         newcol = x - min(x)) # plus or minus? 

state_points <- cbind(st_centroid(states), st_coordinates(st_centroid(states$geometry))) %>% 
  select(X, Y, name, postal) %>% 
  mutate(
    X = ifelse(postal=='NJ', X+0.25, X)) %>% 
  filter(!postal == "DC") # nudge for plotting 

ggmap <- ggplot(data = states) +
  geom_rect(data=recs, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="black",  alpha=0.5) +
  geom_sf(fill = "grey92") +
  # geom_point(data=hauldat, aes(x=lon, y=lat), color="black", fill="black", size=0.1) +
  xlab("Longitude") + ylab("Latitude") +
  annotation_scale(location = "bl", width_hint = 0.5, pad_x = unit(1, "in")) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(2.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_text(data= state_points,aes(x=X, y=Y, label=postal),
            color = "black", fontface = "bold", check_overlap = FALSE, position="dodge") +
  # geom_rect(data=recs, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), color="black",  alpha=0.5) +
  #  geom_tile(data = recs2, aes(x=meanx, y=y, width=2, height=1, fill=newcol), color="black") +
  scale_fill_gradient(low="darkblue", high="white") +
  coord_sf(xlim = c(-80, -65), ylim = c(34, 45.5), expand = FALSE) +
  annotate(geom = "text", x = -69.5, y = 38, label = "Atlantic Ocean", 
           fontface = "italic", color = "grey22", size = 6)  +
  NULL

ggsave(ggmap, filename=here("results","area_map.png"), width=110, height=110, dpi=600, units="mm")

############
# make hovmoller plot 
############

# calculate isotherms 

isos <- seq(9, 15, 2)  # these are the isotherms we want to plot 

# workflow is to fit models to temperature ~ latitude in every year and then use those to predict latitude given temperature 
# these are intentionally not pooled at all because we are just using them to calculate the location of isotherms in a given year (so we don't want to borrow information across years)

if(run_in_parallel==TRUE) {
  library(foreach)
  library(doParallel)
  
  cl <- makeCluster(n_cores) 
  registerDoParallel(cl)
  
  dat_isos <- foreach(i = c(years, years_proj), .combine = bind_rows, .packages = c("brms","dplyr")) %dopar% 
    
    {
      
      dat_bind <- bind_rows(dat_catchonly, dat_test_catchonly) |> 
        select(btemp, year, lat) |> 
        filter(year==i) 
      
      brm_tmp <- brm(btemp ~ lat, data = dat_bind, family = gaussian(), chains = 4, iter = 2000) # this isn't a mistake, there is a reason the model predicts temperature not latitude. if we do lat ~ temp, the isotherms get very flattened, because the model treats the predictor variable as "known" (and there isn't as much variation in latitude as there is in temperature). if we plot temperature vs latitude, there's a lot more variation in temperature. putting btemp as the x-variable, even though it's intuitive because we want to "plug in" isotherm values later, basically throws all that out because the regression is only trying to deal with variation in y given a fixed x. flipping this, so that we're predicting temperature with latitude, is more statistically correct (we're predicting the unknown, variable thing with the known, fixed thing) and lets the model accurately estimate the central tendency and spread of temperatures at a given latitude. then we just invert the formula for a linear regression below to pull out isotherm latitudes. 
      brm_tmp_draws <- as_draws_df(brm_tmp)
      
      tibble(
        isotherm = isos,
        lat = vapply(isos, function(j) {
          median((j - brm_tmp_draws$b_Intercept) / brm_tmp_draws$b_lat)
        }, numeric(1)),
        year = i
      )
      
    }
  
  stopCluster(cl)
  
}
write_rds(dat_isos, here("processed-data","sbt_isotherms.rds"))


if(run_in_parallel==FALSE) {
  dat_isos <- NULL
  out <- NULL
  for(i in c(years, years_proj)) {
    dat_bind <- bind_rows(dat_catchonly, dat_test_catchonly) |> 
      select(btemp, year, lat) |> 
      filter(year==i) 
    
    brm_tmp <- brm(btemp ~ lat, data = dat_bind, family = gaussian(), chains = 4, iter = 2000)
    brm_tmp_draws <- as_draws_df(brm_tmp)
    
    for(j in isos){
      pred <- (j - brm_tmp_draws$b_Intercept) / brm_tmp_draws$b_btemp 
      out$isotherm <- j
      out$lat <- median(pred)
      out$year <- i
      dat_isos <- bind_rows(dat_isos, out)
    }
  }
  
  write_rds(dat_isos, here("processed-data","sbt_isotherms.rds"))
}

# did the isotherms shift north? 
dat_isos %>%
  group_by(isotherm) %>%
  nest() %>%
  mutate(
    model = purrr::map(
      data,
      ~ brm(
        lat ~ year,
        data = .x,
        family = gaussian(),
        chains = 4,
        iter = 2000,
        refresh = 0
      )
    ),
    tidymodel = purrr::map(
      model,
      ~ broom.mixed::tidy(.x, effects = "fixed")
    )
  ) %>%
  unnest(tidymodel) %>%
  dplyr::select(-data, -model)


# get bottom temperature data for plot 
dat_btemp <- bind_rows(dat_catchonly, dat_test_catchonly) %>% 
  mutate(lat_floor_offset = floor(lat) + 0.5) %>% # offset for plotting 
  group_by(year, lat_floor_offset) %>% 
  summarise(btemp = mean(btemp, na.rm=TRUE),
            n = length(haulid[!is.na(btemp)]))

gg_btemp <- ggplot() + 
  geom_raster(data=dat_btemp, aes(x=year, y=lat_floor_offset, fill=btemp)) + 
  scale_fill_gradientn(colors=c("blue3","darkturquoise", "gold", "orangered", "red3"), limits = c(6,24), breaks = c(seq(6, 24, 2)),
                       guide = guide_colourbar(nbin=100, draw.ulim = FALSE, draw.llim = FALSE)) + 
  geom_line(data=dat_isos, aes(x=year, y=lat, group=isotherm)) +
  geom_dl(data=dat_isos, aes(x=year, y=lat, group=isotherm, label=isotherm), method=list(dl.trans(x = x + 0.5, y = y + 0.5, cex=0.8), "first.points")) +
  scale_x_continuous(limits = c(1971, 2017), breaks=seq(1972, 2016, 4), expand = c(0, 0)) +
  scale_y_continuous(limits = c(34.5, 45.5), breaks=seq(35, 45, 1), expand = c(0, 0)) +
  labs(x="Year", y="Latitude") +
  theme_bw() +
  theme(text=element_text(family="sans",size=12,color="black"),
        legend.text = element_text(size=12),
        axis.text=element_text(family="sans",size=8,color="black"), 
        axis.text.x = element_text(angle = 90, hjust = 1), 
        axis.title=element_text(family="sans",size=12,color="black"),
        plot.margin=margin(t = 15, r = 5, b = 5, l = 5, unit = "pt")) + 
  guides(fill = guide_colourbar(barwidth = 0.75, barheight = 25, title="SBT", reverse=TRUE)) +
  NULL

ggsave(gg_btemp, filename=here("results","btemp_lat_time.png"), width=110, height=70, scale = 2.1, dpi=600, units="mm")


############
# make time-series tileplot 
############

#  plot entire time-series
gg_observed <- dat_test_dens %>%
  mutate(Year = (year + min(years_proj) - 1), Latitude = lat_floor, Density=mean_dens, .keep="none") %>%
  bind_rows(dat_train_dens |> mutate(Year = (year + min(years) - 1), Latitude = lat_floor, Density=mean_dens, .keep="none")) |> 
  ggplot() +
  geom_tile(aes(x=Year, y=Latitude, fill=Density)) +
  geom_line(data = time_series_dat, aes(x=Year, y=value, group=feature), color="white") + 
  geom_vline(aes(xintercept = 2006.5), color="white", linetype = "dashed") +
  geom_label(aes(x=Year, y=value, label = feature),
             data = time_series_dat %>% filter(Year == 2002 & !feature=="Warm edge"),
             nudge_y = 0.5,
             size = 4, 
             label.size = 0, 
             fill = "transparent",
             color = "white")  +
  geom_label(aes(x=Year, y=value, label = feature), # moving just this one label downward 
             data = time_series_dat %>% filter(Year == 2002 & feature=="Warm edge"),
             nudge_y = -1.6,
             size = 4, 
             label.size = 0, 
             fill = "transparent",
             color = "white")  +
  scale_x_continuous(breaks=seq(min(years), max(years_proj), 4), limits=c(1971.49, 2016.51)) +
  scale_y_continuous(breaks=seq(min(patches), max(patches), 1), limits=c(34.49, 44.51)) +
  scale_fill_viridis_c() + 
  guides(fill = guide_colorbar(theme = theme(legend.direction = "horizontal"))) +
  theme(legend.position = c(0.3, 0.83), 
        axis.text.x = element_text(angle = 45, vjust=0.8)) + 
  coord_cartesian(expand = FALSE) +
  NULL
ggsave(gg_observed, filename=paste0(here("results"),"/tileplot_timeseries.png"), dpi=600, units="mm", width=75, height=50, scale = 1.7)


############
# make model comparison (bias vs rmse) plot
############


dat_mod_compare <- dat_forecasts_summ  %>% 
  mutate(feature = case_match(feature, "centroid" ~ "Centroid", "warm_edge" ~ "Warm Edge", "cold_edge" ~ "Cold Edge", .default=feature),
         type = ifelse(str_detect(name, "DRM"),"DRM", name))

# name_order <- dat_mod_compare |> 
#   filter(metric=="Bias") |> 
#   group_by(name) |> 
#   summarise(avg = mean(abs(value))) |> 
#   arrange(avg) |> 
#   pull(name)

gg_mod_compare <- dat_mod_compare %>%  
  mutate(
    metric = factor(metric, levels = c("RMSE","Bias")),
    name = factor(name, levels = c("Persistence", "GAM", "DRM T-movement", "DRM T-mortality","DRM T-recruit","DRM null" )),
    feature = factor(feature, levels = c("Cold Edge", "Centroid", "Warm Edge"))
    #   name = reorder_within(name, value, list(metric, feature)) # this option ranks each subplot by its score so the points are sequential in value along the plot
  ) %>%  
  ggplot( )+
  geom_point(aes(name, y=value, shape = type)) + 
  geom_hline(yintercept=0, linetype="dashed") + 
  coord_flip() + 
  #  scale_x_reordered() +
  scale_y_continuous(expand = c(0,0.4)) +
  theme_bw() + 
  facet_grid(feature ~ metric, scales="free", axes="all_y", axis.labels = "all_y") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45),
        legend.position = "none") +
  ggh4x::facet_nested_wrap(~metric+feature, scales = "free") +
  ggh4x::facetted_pos_scales(
    x = list(COL == 2 ~ scale_x_discrete(guide = 'none'),
             COL == 3 ~ scale_x_discrete(guide = 'none'))
  ) +
  NULL

ggsave(gg_mod_compare, filename=here("results","bias_v_rmse.png"), width=110, height=80, dpi=600, units="mm", scale=1.5)


############
# make best DRM time-series plots  
############
mods_to_show <- c("v0.40", "v0.41")

generate_drm_ribbon_plots <- function(modname){
  results_path <-here("results",paste0(modname))
  drm_dens_proj <- read_rds(file.path(results_path, "density_obs_proj.rds")) %>% 
    mutate(year = year + min(years_proj) - 1,
           patch = patch + min(patches) - 1) %>% 
    filter(year < 2017)
  
  drm_edges <- read_rds(file.path(results_path, "range_quantiles_proj.rds")) %>% 
    filter(range_quantiles_proj < Inf,
           year < 2017) %>% 
    mutate(year = year + min(years_proj) - 1,
           range_quantiles_proj = range_quantiles_proj + min(patches)) 
  
  drm_centroids <- read_rds(file.path(results_path, "centroid_proj.rds"))%>% 
    mutate(year = year + min(years_proj) - 1) %>% 
    filter(year < 2017)
  
  centroids_for_ribbon_plot <- drm_centroids |> 
    group_by(year) |> 
    summarise(value_tmp = median(centroid_proj)) |> 
    mutate(feature = "Centroid", name = "DRM") |> 
    bind_rows(points_for_plot %>% filter(feature=='Centroid')) |> 
    rename("Latitude" = value_tmp, "Year" = year) %>% 
    mutate(name = factor(name, levels=c('DRM','Observed','GAM','Persistence')))
  
  centroid_50 <- drm_centroids |> 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.50)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.50)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "50")
  
  centroid_80 <- drm_centroids |> 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.80)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.80)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "80")
  
  centroid_95 <- drm_centroids |> 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.95)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(centroid_proj), prob = 0.95)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "95")
  
  centroid_ribbons <- bind_rows(centroid_50, centroid_80, centroid_95) |> 
    mutate(ci_level = factor(ci_level, levels=c("95","80","50"))) # for plotting
  
  gg_best_drm_centroid <- ggplot() +
    geom_ribbon(data = centroid_ribbons, aes(x=year, ymin=lower, ymax=upper, fill=ci_level), alpha = 0.7) +
    geom_line(data = centroids_for_ribbon_plot, 
              aes(x=Year, y=Latitude, color=name), lwd = 1) + 
    scale_x_continuous(breaks =seq(2007, 2016, 1), limits=c(2007, 2016)) + 
    scale_y_continuous(breaks = seq(37, 39, 1), labels =  seq(37, 39, 1), limits=c(36.5, 40)) +
    scale_color_manual(values=c("black", "#E20134","#8400CD","#009F81"), name="") + 
    #  scale_fill_manual(values=c("#E20134","#8400CD","#009F81"), name="") + 
    scale_fill_brewer(#name = "Credible Interval", 
      guide = "none") +
    labs(x=NULL, y="Latitude", title = "Centroid") + 
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(), 
          plot.title = element_text(hjust = 0.95, vjust = -26.5),
          plot.margin = unit(c(.1,.1,.1,.1), "cm"))
  gg_best_drm_centroid
  
  cold_edges_for_ribbon_plot <- drm_edges |> 
    filter(quantile == 0.95) |> 
    group_by(year) |> 
    summarise(value_tmp = median(range_quantiles_proj)) |> 
    mutate(feature = "Cold Edge", name = "DRM") |> 
    bind_rows(points_for_plot %>% filter(feature=='Cold Edge')) |> 
    rename("Latitude" = value_tmp, "Year" = year) %>% 
    mutate(name = factor(name, levels=c('DRM','Observed','GAM','Persistence')))
  
  cold_edge_50 <- drm_edges |> 
    filter(quantile == 0.95) |> # get cold edge only 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.50)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.50)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "50")
  
  cold_edge_80 <- drm_edges |> 
    filter(quantile == 0.95) |> # get cold edge only 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.80)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.80)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "80")
  
  cold_edge_95 <- drm_edges |> 
    filter(quantile == 0.95) |> # get cold edge only 
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.95)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.95)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "95")
  
 cold_edge_ribbons <- bind_rows(cold_edge_50, cold_edge_80, cold_edge_95) |> 
    mutate(ci_level = factor(ci_level, levels=c("95","80","50"))) # for plotting
  
  gg_best_drm_cold_edge <- ggplot() +
    geom_ribbon(data = cold_edge_ribbons, aes(x=year, ymin=lower, ymax=upper, fill=ci_level), alpha = 0.7) +
    geom_line(data = cold_edges_for_ribbon_plot, 
              aes(x=Year, y=Latitude, color=name), lwd = 1) + 
    scale_x_continuous(breaks =seq(2007, 2016, 1), limits=c(2007, 2016)) + 
    scale_y_continuous(breaks = seq(39, 44, 1), labels =  seq(39, 44, 1), limits=c(39, 44)) +
    scale_color_manual(values=c("black", "#E20134","#8400CD","#009F81"), name="") + 
    #  scale_fill_manual(values=c("#E20134","#8400CD","#009F81"), name="") + 
    scale_fill_brewer(#name = "Credible Interval", 
      guide = "none") +
    labs(x=NULL, y="Latitude", title = "Cold  Edge") + 
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(), 
          plot.title = element_text(hjust = 0.95, vjust = -26.5),
          plot.margin = unit(c(.1,.1,.1,.1), "cm"))
  gg_best_drm_cold_edge
  
  warm_edges_for_ribbon_plot <- drm_edges |> 
    filter(quantile == 0.05) |> 
    group_by(year) |> 
    summarise(value_tmp = median(range_quantiles_proj)) |> 
    mutate(feature = "Warm Edge", name = "DRM") |> 
    bind_rows(points_for_plot %>% filter(feature=='Warm Edge')) |> 
    rename("Latitude" = value_tmp, "Year" = year) %>% 
    mutate(name = factor(name, levels=c('DRM','Observed','GAM','Persistence')))
  
  
  warm_edge_50 <- drm_edges |> 
    filter(quantile == 0.05) |>  
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.50)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.50)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "50")
  
  warm_edge_80 <- drm_edges |> 
    filter(quantile == 0.05) |>   
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.80)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.80)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "80")
  
  warm_edge_95 <- drm_edges |> 
    filter(quantile == 0.05) |>  
    group_by(year) |>
    summarise(
      lower = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.95)
        h[1]
      },
      upper = {
        h <- coda::HPDinterval(coda::as.mcmc(range_quantiles_proj), prob = 0.95)
        h[2]
      },
      .groups = "drop"
    ) |>
    mutate(ci_level = "95")
  
  warm_edge_ribbons <- bind_rows(warm_edge_50, warm_edge_80, warm_edge_95) |> 
    mutate(ci_level = factor(ci_level, levels=c("95","80","50"))) # for plotting
  
  
  
  
  gg_best_drm_warm_edge <- ggplot() +
    geom_ribbon(data = warm_edge_ribbons, aes(x=year, ymin=lower, ymax=upper, fill=ci_level), alpha = 0.7) +
    geom_line(data = warm_edges_for_ribbon_plot, 
              aes(x=Year, y=Latitude, color=name), lwd = 1) + 
    scale_x_continuous(breaks =seq(2007, 2016, 1), limits=c(2007, 2016)) + 
    scale_y_continuous(breaks = seq(35, 38, 1), labels =  seq(35, 38, 1), 
                       limits=c(34.5, 38.5)) +
    scale_color_manual(values=c("black", "#E20134","#8400CD","#009F81"), name="") + 
    #  scale_fill_manual(values=c("#E20134","#8400CD","#009F81"), name="") + 
    scale_fill_brewer(#name = "Credible Interval", 
      guide = "none") +
    labs(x="Year", y="Latitude", title = "Warm  Edge") + 
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, vjust=0.8), 
          plot.title = element_text(hjust = 0.95, vjust = -26.5),
          plot.margin = unit(c(.1,.1,.1,.1), "cm")) 
  gg_best_drm_warm_edge
  
  gg_out <- gg_best_drm_cold_edge + gg_best_drm_centroid + gg_best_drm_warm_edge + plot_layout(ncol=1)
  return(gg_out)}

gg_best_drm <- (wrap_elements(panel = grid::textGrob("DRM T-Recruit")) + wrap_elements(panel = grid::textGrob("DRM T-Mortality"))) / (generate_drm_ribbon_plots(mods_to_show[1]) | generate_drm_ribbon_plots(mods_to_show[2]))+
  plot_layout(heights = c(.05,1,1))

gg_best_drm


ggsave(gg_best_drm, filename=here("results", "best_drm_time.svg"), dpi=600, units="mm", width=75, height=130, scale = 2)
