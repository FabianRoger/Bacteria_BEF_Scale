
# Simulate a BEF experiment

# load the relevant packages
library(here)
library(ggplot2)
library(dplyr)
library(pbapply)

# load the functions
source(here("mcomsimr_simulate_MC2_function.R"))


# set a seed
set.seed(54258748)

# set-up model inputs
species <- 3
patches <- 3
timesteps <- 100
temp_fluc = 0.025
dispersal = 0.05
start_abun = 60
extirp_prob = 0.000001

max_r = 0.5
K_max = 150

int_min = 0
int_max = 1
intra = 1

# landscape parameters
# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75),
             y = c(50, 50, 50))

# generate a random dispersal matrix
d.1 <- mcomsimr::dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

l.test <- 
  pblapply(1:1000, function(z) {
  
  # competition matrices
  si.1 <- matrix(runif(n = species*species, min = int_min, max = int_max), 
                 nrow = species, ncol = species)
  si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
  diag(si.1) <- intra
  
  # species traits
  optima = runif(n = species, 0.2, 0.8)
  env_niche_breadth = runif(n = 1, 0.1, 0.4)
  t.1 <- 
    data.frame(species = 1:species,
               optima = optima,
               env_niche_breadth = env_niche_breadth,
               max_r = max_r,
               K_max = K_max)
  
  # environmental heterogeneity
  e.1 <- Simulate_env(patches = patches, timesteps = timesteps,
                      start_env = c(0.3, 0.5, 0.7), temp_fluc = 0.01)
  
  # Dispersal
  
  # simulate monocultures and mixtures in identical environmental conditions
  MC_disp <- sim_metacomm_BEF(species = species, patches = patches,
                              dispersal = dispersal,
                              timesteps = timesteps, 
                              start_abun = start_abun,
                              extirp_prob = extirp_prob,
                              landscape = l.1, 
                              disp_mat = d.1, 
                              env.df = e.1, 
                              env_traits.df = t.1, 
                              int_mat = si.1,
                              meas_error = 5
  )
  
  # extract the correctly processed data
  MC_disp <- 
    MC_disp$BEF_slope %>%
    filter(time == max(time)) %>%
    group_by(mono_mix, SR) %>% 
    summarise(biomass = mean(biomass))
  
  # extract the BEF-slope
  lm.x <- lm(biomass ~ SR, data = MC_disp)
  B1 <- coef(lm.x)[2]
  names(B1) <- NULL
  
  # get transgressive overyielding
  TO1 <- 
    MC_disp %>%
    group_by(SR) %>%
    summarise(biomass = max(biomass)) %>%
    pull(biomass)
  TO1 <- TO1[2]/TO1[1]
  
  
  # No dispersal
  
  # simulate monocultures and mixtures in identical environmental conditions
  MC_no_disp <- sim_metacomm_BEF(species = species, patches = patches,
                                 dispersal = 0,
                                 timesteps = timesteps, 
                                 start_abun = start_abun,
                                 extirp_prob = extirp_prob,
                                 landscape = l.1, 
                                 disp_mat = d.1, 
                                 env.df = e.1, 
                                 env_traits.df = t.1, 
                                 int_mat = si.1,
                                 meas_error = 5
  )
  
  # extract the correctly processed data
  MC_no_disp <- 
    MC_no_disp$BEF_slope %>%
    filter(time == max(time))
  
  # calculate the species specialisation index
  SI_index <- 
    MC_no_disp %>%
    filter(mono_mix != "mixture") %>%
    group_by(patch) %>%
    filter(biomass == max(biomass)) %>%
    pull(mono_mix) %>%
    unique() %>%
    length()
  SI_index <- SI_index/species
  
  B2 <- 
    sapply(split(MC_no_disp, MC_no_disp$patch), function(x) {
      lm.x <- lm(biomass ~ SR, data = x)
      B <- coef(lm.x)[2]
      names(B) <- NULL
      return(B)
    }) %>%
    mean(.)
  
  TO2 <- 
    sapply(split(MC_no_disp, MC_no_disp$patch), function(x) {
      y <- 
        x %>%
        group_by(SR) %>%
        summarise(biomass = max(biomass)) %>%
        pull(biomass)
      y[2]/y[1]
    }) %>%
    mean(.)
  
  return(data.frame(Scale = c("Large", "Small"),
                    Ave_comp = mean(si.1[si.1 != 1]),
                    SI_index = SI_index,
                    BEF_slope = c(B1, B2),
                    TO = c(TO1, TO2) ) )
  
})

bind_rows(l.test, .id = "run") %>%
  ggplot(data = .,
         mapping = aes(x = Scale, y = BEF_slope)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(width = 0.2) +
  theme_classic()

bind_rows(l.test, .id = "run") %>%
  ggplot(data = .,
         mapping = aes(x = Scale, y = TO)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot(width = 0.2) +
  theme_classic()

# does the difference between large and small scales depend on species specialisation
bind_rows(l.test, .id = "run") %>%
  group_by(run) %>%
  summarise(TO_diff = first(TO) - last(TO),
            SI = mean(SI_index),
            Comp = mean(Ave_comp)) %>%
  filter(Comp > 0.5) %>%
  ggplot(data = .,
         mapping = aes(x = SI, y = TO_diff, colour = Comp)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm") +
  theme_classic()
 





 