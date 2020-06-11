
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale (literature synthesis code)

# load libraries
library(tidyverse)
library(viridis)
library(broom)
library(corrplot)
library(here)

rm(list = ls())

# make a folder to export figures
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# create customised plotting theme
theme_meta <- function(base_size = 12, base_family = "") {
  theme(panel.background =  element_rect(fill = "white"), 
        panel.border =      element_rect(fill="NA", color="black", size=0.35, linetype="solid"),
        axis.line.x = element_line(color="black", size = 0.2),
        axis.line.y = element_line(color="black", size = 0.2),
        panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        axis.ticks.length = unit(-0.16, "cm"),
        axis.title.x = element_text(colour ="black", size = 10, face = "plain", margin=margin(0,0,0,0,"pt")),
        axis.title.y = element_text(colour = "black", size = 10, face = "plain", margin=margin(0,5,0,0,"pt")),
        axis.text.x = element_text(colour = "black", size=10, face = "plain",  margin=margin(10,0,0,0,"pt")),
        axis.text.y = element_text(colour ="black", size=10, face = "plain", margin=margin(0,10,0,0,"pt")),
        axis.ticks.x = element_line(colour = "black", size = 0.4),
        axis.ticks.y = element_line(colour = "black", size = 0.4))
}  


### download data

# the raw data are archived on [Figshare](https://doi.org/10.6084/m9.figshare.12287303.v1)

# citation: Gamfeldt, Lars; Roger, Fabian; Palm, Martin; Hagan, James; Warringer, Jonas; Farewell, Anne (2020): Biodiversity and ecosystem functioning across gradients in spatial scale (literature synthesis). figshare. Dataset. https://doi.org/10.6084/m9.figshare.12287303.v1

# load the raw data
meta_dat_raw <- read_csv( url("https://ndownloader.figshare.com/files/22647539") )

# the raw data file has extra columns in the csv file structure which we remove
meta_dat_raw <- 
  meta_dat_raw %>% 
  select(-contains("X3"))


### based on selection criteria, we removed certain data points from our original raw file

# remove 'mixture best' data points in Fridley (2003)
meta_dat_raw <- filter(meta_dat_raw, Mixture_treatment != "mixture_best")

# remove data from Fox (2002) because it is a pure resource manipulation
meta_dat_raw <- filter(meta_dat_raw, Reference != "Fox_2002")

# remove treatment manipulations that relied on disturbance
meta_dat_raw <- filter(meta_dat_raw, Env_type_1 != "disturbance")


### clean data

# Create a unique identifier for each experiment
meta_dat_raw <- mutate(meta_dat_raw, Experiment_ID = paste(Reference, Experiment_number, sep = "_"))

# Translate ecosystem function values to positive if low value indicates high ecosystem function
meta_dat_raw <- 
  meta_dat_raw %>% 
  group_by(Experiment_ID) %>% 
  mutate(ef_min = min(Ecosystem_function_mean)) %>% ungroup() %>%
  mutate(ef_min = if_else(ef_min < 0, (-ef_min), 0)) %>% 
  mutate(Ecosystem_function_mean = (Ecosystem_function_mean + ef_min)) %>%
  select(-ef_min)


### define functions for the analysis

# standardise_which(x, type)

# standardise data in ways that are not readily available from R
# x = vector of values, type = standardisation type: "max", "max_min" or "none"
standardise_which <- function(x, type = "max") {
  if (type == "max_min") {
    ((x)-min(x))/(max(x)-min(x))
  } else if (type == "max") {
    (x/max(x))
  } else if (type == "none") {
    (x)
  }
}

# bef_slopes_func(data, combs) 

# calculate slope between diversity and function with different numbers of environments
# data = data, combs = vector of environment combinations
bef_slopes_func <- function(data, combs) {
  bind_cols(Experiment_ID = unique(pull(data, Experiment_ID)), 
            bef_slope_env = lm(Ecosystem_function_mean ~ Species_number,
                               filter(data, Environment %in% combs))$coefficients[[2]],
            Environment = paste(combs, collapse = ""))
}

# multiple_hab(data, combs, method)

# calculate the monoculture and mixture values for the different landscape combinations for average overyielding or transgressive overyielding
# data = data, combs = vector of environment combinations, method = normal overyielding "overyield" or transgressive overyielding "trans_overyield"
multiple_hab <- function(data, combs, method) {
  data_combs <- filter(data, Environment %in% combs) %>%
    group_by(Experiment_ID, Mixture_treatment, Mixture_ID) %>%
    summarise(Ecosystem_function_mean = sum(Ecosystem_function_mean)) %>% ungroup() %>%
    group_by(Experiment_ID, Mixture_treatment)
  if (method == "trans_overyield") {
    data_combs %>% summarise(max_mixture_treat = max(Ecosystem_function_mean)) %>% ungroup() %>%
      spread(Mixture_treatment, max_mixture_treat) %>%
      mutate(Environment = paste(combs, collapse = ""))
  } else if (method == "overyield") {
    data_combs %>% summarise(mean_mixture_treat = mean(Ecosystem_function_mean)) %>% ungroup() %>%
      spread(Mixture_treatment, mean_mixture_treat) %>%
      mutate(Environment = paste(combs, collapse = ""))
  }
  else {paste("error")}
}

# meta_stand(data, stand_level, stand)

# Function to standardise the data by within-habitat (all weighted equally) or within-experiment
# data = data, stand_level = level of standardisation "experiment" or "habitat", stand = standardisation type "max", "max_min", "none"
meta_stand <- function(data, stand_level, stand = "max_min") {
  if (stand_level == "experiment") {
    data %>% group_by(Experiment_ID) %>%
      mutate(Ecosystem_function_mean = standardise_which(Ecosystem_function_mean, type = stand)) %>%
      ungroup() 
  } else if (stand_level == "habitat") {
    data %>% group_by(Experiment_ID, Environment) %>%
      mutate(Ecosystem_function_mean = standardise_which(Ecosystem_function_mean, type = stand)) %>%
      ungroup()
  } else {paste("error")}
}

# overyield_func(data, func) 

# calculate overyielding based on mixture and monoculture values
# absolute difference ("abs_diff"): mixture - monoculture
# ratio ("ratio"): mixture/monoculture
# percent change ("percent"): (mixture-monoculture/monoculture*100)

overyield_func <- function(data, func) {
  if (func == "abs_diff") {
    data %>% mutate(overyield_met = (mixture-monoculture))
  } else if (func == "ratio") {
    data %>% mutate(overyield_met = (mixture/monoculture))
  } else if (func == "percent") {
    data %>% mutate(overyield_met = (((mixture/monoculture)-1)*100))
  } else {paste("error")}
}


### create a function to do the analysis with different parameters

# meta_scale_anal(dataset, standardisation_level, standardisation_type, oy_metric, oy_func_calc)

# dataset = data
# standardisation_level = "experiment" or "habitat"
# standardisation_type = "max", "max_min" or "none" from the standardise_which() function
# oy_metric = "abs_diff", "ratio", "percent" from overyield_func() function
# oy_func_calc = "overyield", trans_overyield" from multiple_hab()

meta_scale_anal <- function(dataset = meta_dat_raw, 
                            standardisation_level,
                            standardisation_type,
                            oy_metric,
                            oy_func_calc) {
  
  meta_dat <- meta_stand(data = dataset, stand_level = standardisation_level, stand = standardisation_type)
  
  # for each study and environmental combination calculate:
  # (1) the bef slope
  # (2) average overyielding
  # (3) transgressive overyielding
  
  over_yield_dat <- vector("list", length(unique(meta_dat$Experiment_ID)))
  names(over_yield_dat) <- unique(meta_dat$Experiment_ID)
  env_within <- vector("list")
  bef_slope_within <- vector("list")
  
  for (nm in unique(meta_dat$Experiment_ID)) {
    for (i in seq_along(1:length(unique(meta_dat[meta_dat$Experiment_ID == nm,]$Environment)))) {
      
      combins_env <- combn(x = seq_along(1:length(unique(meta_dat[meta_dat$Experiment_ID == nm,]$Environment))), 
                           m = i, simplify = TRUE)
      dat_ref <- lapply(seq_len(ncol(combins_env)), function(z) meta_dat[meta_dat$Experiment_ID == nm,])
      dat_combs <- lapply(seq_len(ncol(combins_env)), function(y) combins_env[,y])
      
      env_within[[i]] <- bind_rows(mapply(multiple_hab, dat_ref, dat_combs, oy_func_calc, SIMPLIFY = FALSE))
      
      bef_slope_within[[i]] <- bind_rows(mapply(bef_slopes_func, dat_ref, dat_combs, SIMPLIFY = FALSE))
      
      over_yield_dat[[nm]] <- full_join(bind_rows(env_within), bind_rows(bef_slope_within),
                                        by = c("Experiment_ID", "Environment"))
    }
  }
  
  # collapse the list into a dataframe
  
  bef_scale_dat <- bind_rows(over_yield_dat, .id = "Group") %>% 
    filter(Group == Experiment_ID) %>%
    mutate(Number_of_environments = nchar(Environment)) %>%
    mutate(monoculture = (monoculture/Number_of_environments), 
           mixture = (mixture/Number_of_environments)) %>%
    select(Experiment_ID, Environment, Number_of_environments, monoculture, mixture, bef_slope_env)
  
  
  # for each experiment, calculate the slope between scale and:
  # (1) the bef slope
  # (2) average overyielding (average or transgressive depending on argument)
  
  slopes_dat <- overyield_func(bef_scale_dat, func = oy_metric)
  experiments <- unique(slopes_dat$Experiment_ID)
  
  bef_scale_slopes <- vector("list", length(experiments))
  for (i in seq_along(experiments)) {
    meta_lm1 <- lm(overyield_met ~ Number_of_environments,
                   data = slopes_dat[slopes_dat$Experiment_ID == experiments[i],])
    out1 <- summary(meta_lm1)
    
    meta_lm2 <- lm(bef_slope_env ~ Number_of_environments,
                   data = slopes_dat[slopes_dat$Experiment_ID == experiments[i],])
    out2 <- summary(meta_lm2)
    
    bef_scale_slopes[[i]] <- tibble(Est = out1$coefficients[c(2)],
                                    se_est = out1$coefficients[c(4)],
                                    Est_bef_slope = out2$coefficients[c(2)])
  }
  
  # organise the output into a list
  
  list(standardisation_exp_hab = standardisation_level,
       overyield_metric = oy_metric,
       normal_trans_oy = oy_func_calc,
       raw_data = slopes_dat,
       slopes_data = bind_cols(Stand = rep(standardisation_level, length(experiments)),
                               Stand_type = rep(standardisation_type, length(experiments)),
                               Metric = rep(oy_metric, length(experiments)),
                               oy_function = rep(oy_func_calc, length(experiments)),
                               Experiment_ID = experiments, Slopes = bind_rows(bef_scale_slopes)$Est, 
                               Slopes_se = bind_rows(bef_scale_slopes)$se_est,
                               Slopes_bef = bind_rows(bef_scale_slopes)$Est_bef_slope))
  
}


# meta_scale_out(data, data_output)

# function that to output data from meta_scale_anal() allowing easily varying parameters
# data = data
# data_output = "slopes_data" or "raw_data"
# run_list = list with four parameters corresponding to standardisation_level, standardisation_type, oy_metric and oy_func_calc

meta_scale_out <- function(data = meta_dat_raw, data_output = "slopes_data", run_list = run_list) {
  meta_scale_runs <- vector("list")
  for (i in seq_along(1:(length(run_list)))) {
    data_out <- meta_scale_anal(dataset = data, 
                                standardisation_level = run_list[[i]][1],
                                standardisation_type = run_list[[i]][2],
                                oy_func_calc = run_list[[i]][3],
                                oy_metric = run_list[[i]][4])
    
    meta_scale_runs[[i]] <- data_out[[data_output]]
  }
  
  bind_rows(meta_scale_runs, .id = "run")
}


### calculate explanatory variables and bind them into a dataframe

# total average overyielding: 'mean_ave_overyield'
# coefficient of variation in average overyielding among environments: 'cv_ave_overyield' 
# coefficient of variation in average ecosystem function across environments: 'habitat_cv'  

hab_oy <- 
  meta_dat_raw %>%
  group_by(Experiment_ID, Environment) %>%
  mutate(habitat_mean = mean(Ecosystem_function_mean)) %>% 
  ungroup() %>% group_by(Experiment_ID) %>%
  mutate(habitat_cv = ((sd(habitat_mean)/mean(habitat_mean))*100)) %>%
  ungroup() %>%
  group_by(Experiment_ID, Environment, Mixture_treatment) %>%
  summarise(mixture_treatment_mean = mean(Ecosystem_function_mean),
            habitat_mean = mean(habitat_mean),
            habitat_cv = mean(habitat_cv)) %>%
  spread(key = Mixture_treatment, value = mixture_treatment_mean) %>%
  mutate(ave_overyield = mixture/monoculture) %>% 
  ungroup() %>%
  group_by(Experiment_ID) %>%
  mutate(Pearson_r = if_else(cor(habitat_mean, ave_overyield) > 0, 1, -1)) %>%
  ungroup() %>%
  group_by(Experiment_ID) %>%
  summarise(habitat_mean = mean(habitat_mean),
            habitat_cv = mean(habitat_cv),
            mean_ave_overyield = mean(ave_overyield),
            cv_ave_overyield = ((sd(ave_overyield)/mean(ave_overyield))*100),
            Pearson_r_hab_oy = mean(Pearson_r)) %>%
  gather(key = "hab_var", value = "hab_value",
         habitat_cv, cv_ave_overyield) %>%
  group_by(Experiment_ID) %>%
  mutate(hab_oy_metric = (mean(hab_value)-sd(hab_value))) %>%
  ungroup() %>%
  spread(key = hab_var, value = hab_value) %>%
  mutate(hab_oy_metric = hab_oy_metric*Pearson_r_hab_oy)


### calculate species specialisation index for each experiment

ss_dat <- 
  meta_dat_raw %>% 
  filter(Mixture_treatment != "mixture") %>% 
  group_by(Experiment_ID, Environment) %>%
  mutate(max_eco_func = max(Ecosystem_function_mean)) %>% 
  ungroup() %>%
  mutate(n_spp = if_else(Ecosystem_function_mean == max_eco_func, Mixture_ID, "none")) %>%
  filter(n_spp != "none") %>% 
  group_by(Experiment_ID) %>%
  summarise(n_spp_best = length(unique(n_spp)),
            n_envs = length(unique(Environment))) %>% 
  ungroup() %>%
  mutate(n_spp_env_prop = n_spp_best/n_envs) %>% 
  mutate(species_sorting = (n_spp_env_prop - (1/n_envs))/(1-(1/n_envs)))


# bind these explanatory variables into a dataframe with moderators

explan_vars <- 
  full_join(meta_dat_raw %>%
  select(Experiment_ID, Realm, Trophic_level, Lineage, Env_type_1, Env_type_2, 
         Ecosystem_function_type, Environment) %>%
  group_by(Experiment_ID) %>%
  mutate(Environment = max(Environment)) %>% ungroup() %>%
  distinct() %>%
  select(-Environment),
  full_join(hab_oy, ss_dat, by = "Experiment_ID"), 
  by = "Experiment_ID")


### run the meta_scale_out() function and join this data to the explanatory variables

# set up a list of runs for this function

# these are the arguments that will be called from meta_scale_anal() function:

# standardisation_level = "experiment" or "habitat"
# standardisation_type = "max", "max_min" or "none" from the standardise_which() function
# oy_metric = "abs_diff", "ratio", "percent" from overyield_func() function
# oy_func_calc = "overyield", trans_overyield" from multiple_hab() function

# write these arguments into a list
run_list <- list(run_1 = c("experiment", "max", "trans_overyield", "ratio"),
                 run_2 = c("experiment", "max", "overyield", "ratio"))


# run the meta_scale_out() function to output raw values for each experiment
meta_scale_raw <- 
  full_join(meta_scale_out(data = meta_dat_raw, data_output = "raw_data", run_list = run_list), 
                            explan_vars, by = "Experiment_ID")

# clean this output data
meta_scale_raw <- 
  meta_scale_raw %>% 
  select(-monoculture, -mixture) %>%
  spread(key = run, value = overyield_met) %>%
  rename(trans_overyield = "1", average_overyield = "2", bef_slope = bef_slope_env) %>%
  gather(key = effect_type, value = effect_value, 
         bef_slope, average_overyield, trans_overyield)


# check the meta_scale_raw calculations manually for different experiments

# select the experiment id
exp_test <- c("Gamfeldt_et_al_2005_one")

# select the environments
env_test <- c(1, 2)
env_id <- paste(env_test, collapse = "")

# bef-slope for each combination of environments

lm(Ecosystem_function_mean ~ Species_number, 
   data = meta_dat_raw %>% 
     filter(Experiment_ID == exp_test) %>%
     mutate(Ecosystem_function_mean = standardise_which(Ecosystem_function_mean)) %>% 
     filter(Environment %in% env_test))

# compare with output from the function

meta_scale_raw %>% 
  filter(Experiment_ID == exp_test) %>%
  select(Environment, Number_of_environments, effect_type, effect_value) %>%
  filter(Environment == env_id, 
         effect_type == "bef_slope")

# average overyielding for each combination of environments
meta_dat_raw %>% 
  filter(Experiment_ID == exp_test) %>%
  mutate(Ecosystem_function_mean = standardise_which(Ecosystem_function_mean)) %>%
  filter(Environment %in% env_test) %>%
  group_by(Species_number, Mixture_treatment, Mixture_ID) %>%
  summarise(Ecosystem_function_mean = sum(Ecosystem_function_mean)) %>%
  ungroup() %>%
  group_by(Species_number, Mixture_treatment) %>%
  summarise(max_ef = mean(Ecosystem_function_mean)) %>%
  ungroup() %>%
  select(-Species_number) %>%
  spread(Mixture_treatment, max_ef) %>%
  mutate(oy_perc = (((mixture/monoculture)-1)*100),
         oy_ratio = mixture/monoculture)

# compare with output from the function

meta_scale_raw %>% 
  filter(Experiment_ID == exp_test) %>%
  select(Environment, Number_of_environments, effect_type, effect_value) %>%
  filter(Environment == env_id, 
         effect_type == "average_overyield")

# transgressive overyielding for each combination of environments

meta_dat_raw %>% 
  filter(Experiment_ID == exp_test) %>%
  mutate(Ecosystem_function_mean = standardise_which(Ecosystem_function_mean)) %>%
  filter(Environment %in% env_test) %>%
  group_by(Species_number, Mixture_treatment, Mixture_ID) %>%
  summarise(Ecosystem_function_mean = sum(Ecosystem_function_mean)) %>%
  ungroup() %>%
  group_by(Species_number, Mixture_treatment) %>%
  summarise(max_ef = max(Ecosystem_function_mean)) %>%
  ungroup() %>%
  select(-Species_number) %>%
  spread(Mixture_treatment, max_ef) %>%
  mutate(to_perc = (((mixture/monoculture)-1)*100),
         to_ratio = mixture/monoculture)

# compare with output from the function

meta_scale_raw %>% 
  filter(Experiment_ID == exp_test) %>%
  select(Environment, Number_of_environments, effect_type, effect_value) %>%
  filter(Environment == env_id, 
         effect_type == "trans_overyield")

# checked enough that we are confident? yes.


### table S1

# this table lists all studies used in the literature synthesis and some additional metadata
# slight aesthetic changes were made for the manuscript

# Set up unique names for the experimental id's
match_names <- tibble(Experiment_ID = sort(unique(meta_scale_raw$Experiment_ID)),
                      Experiment_code = LETTERS[1:length(unique(meta_scale_raw$Experiment_ID))])

table_s1 <- left_join(match_names, 
                      select(meta_scale_raw, Experiment_ID, Realm, n_envs, Trophic_level, Lineage),
                      by = "Experiment_ID") %>%
  distinct() %>% 
  select(Experiment_code, Experiment_ID, Realm, n_envs, Trophic_level, Lineage)

write_csv(table_s1, here("figures/table_s1.csv") )


### fig. S7, 8, 9

# in these graphs, we plot the relationship between scale and: 
# (1) the bef-slope (fig. S7); 
# (2) average overyielding (fig. S8) and; 
# (3) transgressive overyielding (fig. S9) for each study included in the literature synthesis

effects <- unique(meta_scale_raw$effect_type)
effect_names <- c("BEF-slope", "ave-OY", "trans-OY")
figs7_9_names <- c("S7", "S8", "S9")

raw_plots <- vector("list", length = length(effects))

for (i in seq_along(1:length(effects))) {
  raw_plots[[i]] <- ggplot(data = full_join(meta_scale_raw %>% 
                                              filter(effect_type == effects[i]),
                                            match_names, by = "Experiment_ID"),
                           aes(x = Number_of_environments, y = effect_value,
                               colour = Realm, shape = Ecosystem_function_type)) +
    geom_smooth(se = FALSE, method = "lm", colour = "black", size = 1) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_x_continuous(breaks = c(1:8)) +
    facet_wrap(~Experiment_code, scales = "free", ncol = 5) +
    scale_colour_viridis_d() +
    ylab(effect_names[i]) +
    xlab("scale") +
    theme_meta() +
    labs(shape = "Ecosystem function", colour = "Realm") +
    theme(legend.position = c(0.8, 0.075),
          legend.box = "horizontal",
          legend.key = element_rect(fill = "white"))
  
  ggsave(raw_plots[[i]],
         filename = paste( here("figures"),paste(paste("Fig_", figs7_9_names[i], sep = "_"),".png", sep = ""), sep = "/" ),
         dpi = 400, units = "cm",
         width = 19,
         height = 24)
}


# fig. 3Ai-Aiii inset plots

# this plots the relationship between scale and:
# (1) the bef-slope
# (2) average overyielding
# (3) transgressive overyielding
# for a single experiment ("Blake_Duffy_2010_one") as an example

# note that this graph was aesthetically modified for the main manuscript

slope_type <- unique(meta_scale_raw$effect_type)
slope_names <- c("BEF-slope", "ave-OY", "trans-OY")
ref <- c("Blake_and_Duffy_2010_one")
fig3_names <- c("3Ai", "3Aii", "3Aiii")

raw_plot_bd <- vector("list", length(slope_type))
for (i in seq_along(1:length(slope_type))) {
  raw_plot_bd[[i]] <- ggplot(data = meta_scale_raw %>% 
                               filter(Experiment_ID == ref) %>%
                               filter(effect_type == slope_type[i]),
                             aes(x = Number_of_environments, y = effect_value)) +
    geom_smooth(se = FALSE, method = "lm", colour = "black", size = 0.5) +
    geom_jitter(size = 1.5, alpha = 0.8, width = 0.05) +
    scale_x_continuous(breaks = c(1:4)) +
    scale_colour_viridis_d() +
    ylab(slope_names[i]) +
    xlab("scale") +
    theme_meta() +
    theme(axis.title.x = element_text(size = 8, margin = margin(t = 3.5, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 8, margin = margin(t = 0, r = 3.5, b = 0, l = 0)),
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))
  
  ggsave(raw_plot_bd[[i]],
         filename = paste( here("figures"),paste(paste("Fig_", fig3_names[i], sep = "_"),".svg", sep = ""), sep = "/" ),
         units = "cm",
         width = 4.5,
         height = 3.5)
}



# run the meta_scale_out() function to output raw values for each experiment

meta_scale_slope <- full_join(meta_scale_out(data = meta_dat_raw, 
                                             data_output = "slopes_data",
                                             run_list = run_list), 
                              explan_vars, by = "Experiment_ID")

# clean the data to test hypotheses about the distribution of slopes among experiments

meta_scale_slope <- meta_scale_slope %>% 
  select(-oy_function, -Slopes_se) %>%
  spread(key = run, value = Slopes) %>%
  rename(bef_slope = Slopes_bef, trans_oy_slope = "1", average_oy_slope = "2") %>%
  gather(key = slope_type, value = slope_value, bef_slope, trans_oy_slope, average_oy_slope)


### main text statistics

# here, we test whether hypotheses about the distribution of the slopes between scale and:
# (1) the bef slope
# (2) average overyielding
# (3) transgressive overyielding

# for this, we use Wilcoxon signed rank tests

# define the three different slopes
slope_t <- meta_scale_slope$slope_type %>% 
  unique()

# define the alternative hypotheses for each slope
h0_wil <- c("two.sided", "greater", "two.sided")

wt_out <- vector("list", length(slope_t))

for (i in seq_along(1:length(slope_t))) {
  w <- meta_scale_slope %>% 
    filter(slope_type == slope_t[i]) %>%
    pull(slope_value)
  
  wt <- wilcox.test(w, alternative = h0_wil[i], mu = 0, 
                    paired = FALSE, exact = FALSE,
                    conf.int = FALSE)
  
  wt_out[[i]] <- tibble(slope = slope_t[i],
                        test = h0_wil[i],
                        w_stat = wt$statistic,
                        p_val = wt$p.value)
}

wt_out <- bind_rows(wt_out)

# adjust the P-values and output the results into a table
wt_out %>% 
  mutate(p_val_hom = p.adjust(wt_out$p_val, method = "hommel"),
         p_val_BH = p.adjust(wt_out$p_val, method = "BH")) %>% 
  write_csv( here("figures/Main_text_wilcox_tests.csv") )

# manually check these results for consistency
meta_scale_slope %>% 
  filter(slope_type == slope_t[2]) %>%
  pull(slope_value) %>%
  wilcox.test(alternative = h0_wil[3], mu = 0, 
              paired = FALSE, exact = FALSE,
              conf.int = FALSE)


# fig. 3A

# This graph compares the slope between scale and:
# (1) the bef slope
# (2) average overyielding
# (3) transgressive overyielding
# for all experiments included in the literature synthesis

# note that this graph was aesthetically modified for the main manuscript

fig_3A <- ggplot(data = meta_scale_slope %>% 
               mutate(slope_type = factor(slope_type, 
                                          levels(factor(slope_type))[c(2,1,3)])) %>%
               mutate(slope_type = factor(slope_type, 
                                          labels = c("BEF-scale", "ave-OY-scale", "trans-OY-scale"))),
             aes(x = slope_type, y = slope_value)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.35, alpha = 0.6, outlier.shape = NA) +
  geom_jitter(size = 2, alpha = 0.3, width = 0.1, shape = 16) +
  scale_colour_viridis_d(begin = 0, end = 0.9) +
  xlab("") +
  ylab("slope estimate") +
  theme_meta() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(fig_3A,
       filename = here("figures/Fig_3A.svg"),
       units = "cm",
       width = 8,
       height = 8)


### Supplementary material

# table S2

# here, we test whether hypotheses about the distribution of the slopes between scale and:
# (1) the bef slope
# (2) average overyielding
# (3) transgressive overyielding

# using parametric one-sample t-tests as opposed to non-parametric Wilcoxon signed rank tests
# we also check the influence of outliers in the data

# for this, we write a function called t_test_slopes with arguments:
# data = the meta_scale_slope dataset
# slope_type = vector of slope types
# hypotheses = two-sided or one-sided

t_test_slopes <- function(data, slopes, hypotheses) {
  t_tests <- vector("list", length(slopes))
  for (i in seq_along(1:length(slopes))) {
    t_dat <- data %>% 
      filter(slope_type == slopes[i]) %>%
      pull(slope_value)
    t_obj <- t.test(t_dat, mu = 0, alternative = hypotheses[i])
    t_tests[[i]] <- tibble(slope = slopes[i],
                           estimate = t_obj$estimate,
                           se = t_obj$stderr,
                           t_stat = t_obj$statistic,
                           p_val = t_obj$p.value)
  }
  t_tests <- bind_rows(t_tests)
  mutate(t_tests, p_val_hom = p.adjust(t_tests$p_val, method = "hommel"),
         p_val_BH = p.adjust(t_tests$p_val, method = "BH"))
}

# apply this function to the full dataset and write the output into a table
t_test_full <- t_test_slopes(data = meta_scale_slope,
              slopes = c(unique(meta_scale_slope$slope_type)),
              hypotheses = c("two.sided", "two.sided", "greater")) %>%
  mutate(slope = factor(slope, levels = c("bef_slope", "average_oy_slope", "trans_oy_slope"))) %>%
  arrange(slope)

write_csv(t_test_full, here("figures/Table_S2_all_data.csv") )

# check these results for consistency
t.test(filter(meta_scale_slope, slope_type == "bef_slope") %>% 
         pull(slope_value), mu = 0, alternative = c("two.sided"))


# calculate extreme outliers (> 3 X IQR from 0.25 or 0.75 quantiles) for each slope:
# (1) the bef slope
# (2) average overyielding
# (3) transgressive overyielding

outlier_refs <- meta_scale_slope %>% group_by(slope_type) %>%
  mutate(Q1 = quantile(slope_value, probs = c(0.25))[[1]],
         Q3 = quantile(slope_value, probs = c(0.75))[[1]],
         Inter_Q = IQR(slope_value)) %>%
  mutate(lower_fence = Q1 - (3*Inter_Q),
         upper_fence = Q3 + (3*Inter_Q)) %>%
  mutate(outlier = if_else(slope_value > upper_fence, "out",
                           if_else(slope_value < lower_fence, "out", "normal"))) %>%
  ungroup() %>%
  filter(outlier == "out") %>% 
  select(Experiment_ID, slope_type, outlier) %>%
  distinct() %>%
  pull(Experiment_ID) %>% 
  unique()

# get all combinations of these diferent outliers
outlier_ref_comb <- do.call(c,
        lapply(seq_along(outlier_refs), combn, x = outlier_refs , simplify = FALSE)
        )

# repeat the one-sample t-tests but excluding all combinations of these outliers
t_test_outliers <- vector("list", length(outlier_ref_comb ))
for (i in seq_along(1:length(outlier_ref_comb ))) {
  t_test_outliers[[i]] <- t_test_slopes(data = meta_scale_slope %>% 
                  filter(!Experiment_ID %in% outlier_ref_comb[[i]]),
                slopes = unique(meta_scale_slope$slope_type),
                hypotheses = c("two.sided", "two.sided", "greater"))
}

t_test_outliers <- bind_rows(t_test_outliers, .id = "reference_removed")

# get references that were removed
ref_vec <- vector()
for (i in seq_along(1:length(outlier_ref_comb))) {
  ref_vec[i] <- paste(c(outlier_ref_comb[[i]]), collapse = "-")
}

# correct the p-values and output the results
t_test_outliers <- t_test_outliers %>%
  mutate(reference_removed_id = rep(ref_vec, each = 3)) %>%
  mutate(slope = factor(slope, 
                        levels = c("bef_slope", "average_oy_slope", "trans_oy_slope"))) %>%
  arrange(reference_removed, slope)

write_csv(t_test_outliers, here("figures/Table_S2_outliers_removed.csv")  )


### Main text

# we used model selection to explore how four predictor variables: 
# (1) species specialisation index: species_sorting
# (2) total average overyielding: mean_ave_overyield
# (3) coefficient of variation among environments in average overyielding: cv_ave_overyield
# (4) coefficient of variation among environments in average ecosystem functioning: habitat_cv

# affected the slope between scale and transgressive overyielding between experiments

## make sure the data is clean

# check the distribution of these variables
filter(meta_scale_slope, slope_type == "trans_oy_slope") %>%
  select(slope_value, species_sorting, cv_ave_overyield, mean_ave_overyield, habitat_cv) %>%
  gather(key = "variable", value = "value") %>%
  ggplot(aes(x = value, fill = variable)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free") +
  theme_classic()

# check the correlation among these variables
pairs(select(filter(meta_scale_slope, slope_type == "trans_oy_slope"),
             slope_value, species_sorting, cv_ave_overyield, mean_ave_overyield, habitat_cv))

corrplot(cor(select(filter(meta_scale_slope, slope_type == "trans_oy_slope"),
                slope_value, species_sorting, cv_ave_overyield, mean_ave_overyield, habitat_cv)),
         method = "number", type = "lower")


# set up a function to run different models that can then be compared
lm_bef_scale <- function(data = meta_scale_slope, slope = "trans_oy_slope",
                         explans = exp_vars, outliers = NA) {
  cof_out <- vector("list", length(explans))
  names(cof_out) <- seq_along(1:length(explans))
  diag_out <- vector("list", length(explans))
  names(diag_out) <- seq_along(1:length(explans))
  for (i in seq_along(1:length(explans))) {
    lm_dat <- data %>% filter(slope_type == slope) %>%
      filter(!Experiment_ID %in% outliers)
    lm_exp <- lm(reformulate(explans[[i]], "slope_value"),
                 data = lm_dat)
    cof_out[[i]] <- tidy(lm_exp)
    diag_out[[i]] <- glance(lm_exp)
  }
  
  full_join(bind_rows(cof_out, .id = "model"), bind_rows(diag_out, .id = "model"),
            by = "model")
}

# set up the explanatory variables for the different models
exp_vars <- list(c("species_sorting", "mean_ave_overyield", "habitat_cv", "cv_ave_overyield"),
                 c("species_sorting*mean_ave_overyield"), 
                 c("species_sorting*habitat_cv"),
                 c("species_sorting*cv_ave_overyield"),
                 c("species_sorting"),
                 c("cv_ave_overyield"),
                 c("habitat_cv"),
                 c("mean_ave_overyield"))

# run this set of models with all the data
trans_y_lm <- lm_bef_scale(data = meta_scale_slope, 
                           slope = "trans_oy_slope", 
                           explans = exp_vars, outliers = NA)

# arrange the models by AIC and examine the results
View(trans_y_lm %>% arrange(AIC))


# table S3

# in this table, we output a summary of the models fit to the slopes between scale and transgressive overyielding with all the data

# clean this output, extract relevant columns and output the table
trans_y_lm <- 
  trans_y_lm %>% 
  select(model, term, r.squared, AIC) %>%
  filter(term != "(Intercept)") %>%
  group_by(model, r.squared, AIC) %>%
  summarise(terms = paste(term, collapse = "+")) %>% 
  ungroup() %>%
  arrange(AIC) %>% select(terms, r.squared, AIC) %>%
  mutate(delta_AIC =  AIC - (min(AIC))) %>%
  mutate(AIC_wt_start = exp(1)^(-0.5*delta_AIC)) %>%
  mutate(AIC_wt = AIC_wt_start/sum(AIC_wt_start)) %>%
  select(-AIC_wt_start)

write_csv(trans_y_lm, here("figures/Table_S3.csv") )


# run this set of models without the major outlier: Dzialowski and Smith 2008
trans_y_lm_out <- lm_bef_scale(data = meta_scale_slope, 
                               slope = "trans_oy_slope", 
                               explans = exp_vars, outliers = c("Dzialowski_Smith_2008_one"))

# examine the results
View(trans_y_lm_out %>% arrange(AIC))


# table S4

# In this table, we output a summary of the models fit to the slopes between scale and transgressive overyielding with all the data when a major outlier is removed

# clean this output, extract relevant columns and output the table
trans_y_lm_out <- 
  trans_y_lm_out %>% 
  select(model, term, r.squared, AIC) %>%
  filter(term != "(Intercept)") %>%
  group_by(model, r.squared, AIC) %>%
  summarise(terms = paste(term, collapse = "+")) %>% 
  ungroup() %>%
  arrange(AIC) %>% select(terms, r.squared, AIC) %>%
  mutate(delta_AIC =  AIC - (min(AIC))) %>%
  mutate(AIC_wt_start = exp(1)^(-0.5*delta_AIC)) %>%
  mutate(AIC_wt = AIC_wt_start/sum(AIC_wt_start)) %>%
  select(-AIC_wt_start)

write_csv(trans_y_lm_out, here("figures/Table_S4.csv"))


# figure 3B

# from this model selection, the best model by a significant margin is the species sorting*mean model
# we plot this model for figure 3B

# the figure in the publication is modified for aesthetics

# Using the full dataset
ss_lm <- lm(slope_value ~ species_sorting*mean_ave_overyield, 
            data = filter(meta_scale_slope, slope_type == "trans_oy_slope"))

# check model assumptions
plot(ss_lm)

# check model summary
summary(ss_lm)

# refit the model without the major outlier: Dzialowski_Smith_2008_one
ss_lm_out <- lm(slope_value ~ species_sorting*mean_ave_overyield, 
                data = filter(meta_scale_slope, slope_type == "trans_oy_slope",
                              Experiment_ID != "Dzialowski_Smith_2008_one"))

# check the model assumptions
plot(ss_lm_out)

# check model summary
summary(ss_lm_out)


# calculate species_sorting slope at different mean_ave_overyield values

# check the range of mean_ave_overyield values
av_oy_vector <- 
  filter(meta_scale_slope, slope_type == "trans_oy_slope") %>% 
  pull(mean_ave_overyield)
av_oy_vector

# calculate various summary statistics
range(av_oy_vector)
mean(av_oy_vector)
sd(av_oy_vector)

# set values of explanatory variables to derive model predictions for
new_d <- expand.grid(species_sorting = seq(from = 0, to = 1, by = 0.1),
                     mean_ave_overyield = c((mean(av_oy_vector)-sd(av_oy_vector)), 
                                            (mean(av_oy_vector)),
                                            (mean(av_oy_vector)+sd(av_oy_vector))))

# replicate these values to use with the model without the outlier
new_d_out <- new_d

# get model predictions from these values
new_d$slope_value <- predict(ss_lm, new_d)
new_d_out$slope_value <- predict(ss_lm_out, new_d_out)

# Plot the graph with predictions
fig_3B <- ggplot(data = meta_scale_slope %>% 
                   filter(!slope_type %in% c("bef_slope","average_oy_slope")),
       aes(x = species_sorting, y = slope_value, colour = mean_ave_overyield)) +
  geom_hline(yintercept = 0, alpha = 0.2) +
  geom_line(data = new_d, aes(group = (mean_ave_overyield))) +
  geom_line(data = new_d_out, aes(group = (mean_ave_overyield)), 
            linetype = "dashed", alpha = 0.5) +
  geom_jitter(alpha = 1, width = 0.05, shape = 16, size = 2) +
  scale_colour_viridis_c() +
  xlab("species specialisation index") +
  ylab("trans-OY-scale slope") +
  labs(colour = "AO mean") +
  theme_meta() +
  theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5),
        legend.key.size = unit(0.25, "cm"))

ggsave(fig_3B,
       filename = here("figures/Fig_3B.svg"),
       units = "cm",
       width = 12,
       height = 8)



# literature synthesis data assessment

# here, we examine how many studies, which were excluded and why etc.


### download data

# the raw data are archived on [Figshare](https://doi.org/10.6084/m9.figshare.12287303.v1)

# citation: Gamfeldt, Lars; Roger, Fabian; Palm, Martin; Hagan, James; Warringer, Jonas; Farewell, Anne (2020): Biodiversity and ecosystem functioning across gradients in spatial scale (literature synthesis). figshare. Dataset. https://doi.org/10.6084/m9.figshare.12287303.v1

# load the raw data
data_ass <- read_csv( url("https://ndownloader.figshare.com/files/22647536") )

# check the data
View(data_ass)

colnames(data_ass)

# reorder the columns
data_ass <- 
  data_ass %>% 
  select(meta_analysis_database:reason, data_available_y_n:Note_2, inclusion_exclusion)

# check this database
summary(data_ass)

# check for unique values in each column
lapply(data_ass, function(x) unique(x)) 

lapply(data_ass, function(x) length(unique(x))) 

# how many unique publications did we assess for each of the published meta-analyses?
data_ass %>% 
  mutate(unique_id = paste(reference_id, year, journal, sep = ".")) %>%
  group_by(meta_analysis_database) %>%
  summarise(n = n()) %>% ungroup() %>%
  mutate(total_pubs = sum(n))

# how many unique publications did we assess across the published meta-analyses
data_ass %>% 
  mutate(unique_id = paste(reference_id, year, journal, sep = ".")) %>%
  pull(reference_id) %>% 
  unique() %>% 
  length()

# extract studies that were deemed to have a suitable design
filter(data_ass, suitable_design_y_n == "y") %>%
  pull(reference_id) %>% 
  unique() %>% 
  sort()

filter(data_ass, suitable_design_y_n == "y") %>%
  pull(reference_id) %>% 
  unique() %>% 
  sort() %>%
  length()

# extract the studies with a suitable design and for which data were available
filter(data_ass, suitable_design_y_n == "y", data_available_y_n == "y") %>%
  pull(reference_id) %>% 
  unique() %>% 
  sort()

filter(data_ass, suitable_design_y_n == "y", data_available_y_n == "y") %>%
  pull(reference_id) %>% 
  unique() %>% 
  length()

# extract data for studies with a suitable design but for which data were not available
contacts <- 
  filter(data_ass, suitable_design_y_n == "y", data_available_y_n == "n") %>%
  pull(reference_id) %>% 
  unique()
contacts

# Examine which studies were and weren't included from these for which data were not available
filter(data_ass, reference_id %in% contacts) %>%
  select(reference_id, data_available_y_n, contact_authors_y_n, authors_emailed_y_n,
         Note_1, Note_2, inclusion_exclusion) %>%
  View()

# Smith_Allcock_1985: could not find the author's contact details
# Fridley_2002: incorrect data provided after contacting the author
# Nicklaus et al. (2006): authors contacted but did not respond

# De Boeck et al. 2008: author provided the data

# Extract the data for the studies that were actually included
filter(data_ass, inclusion_exclusion == "inclusion") %>%
  pull(reference_id) %>% 
  unique() %>% 
  sort()

filter(data_ass, inclusion_exclusion == "inclusion") %>%
  pull(reference_id) %>% 
  unique() %>% 
  sort() %>% 
  length()

# 25 unique studies overall (26 experiments in total)
