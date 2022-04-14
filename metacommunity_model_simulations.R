
# Simulate a BEF experiment

# maybe I'm getting a bit off track but I need to explore why get different results
# from single patch model compared to considering single patches in a metacommunity
# model with zero dispersal... It doesn't really make sense now

# maybe it's the way the model is implemented that means the uncertainty propagates differently?

# get landscapes and shit

# set a seed
set.seed(54258748)

# set-up model inputs
species <- 3
patches <- 3
timesteps <- 500
temp_fluc = 0.025
dispersal = 0
start_abun = 60
extirp_prob = 0.000001

optima = seq(0.2, 0.8, length.out = species)
env_niche_breadth = c(0.2, 0.2, 0.2)
max_r = 0.5
K_max = 150

int_min = 1
int_max = 1
intra = 1

# landscape parameters
# generate a random landscape
l.1 <- 
  data.frame(x = c(25, 50, 75),
             y = c(50, 50, 50))

# generate a random dispersal matrix
d.1 <- mcomsimr::dispersal_matrix(landscape = l.1, torus = TRUE, kernel_exp = 0.1, plot = FALSE)

# generate species environmental optima
t.1 <- 
  data.frame(species = 1:species,
             optima = optima,
             env_niche_breadth = env_niche_breadth,
             max_r = max_r,
             K_max = K_max)

# competition matrices
si.1 <- matrix(runif(n = species*species, min = int_min, max = int_max), 
               nrow = species, ncol = species)
si.1[lower.tri(si.1)] = t(si.1)[lower.tri(si.1)]
diag(si.1) <- intra
head(si.1)

# environmental heterogeneity
e.1 <- Simulate_env(patches = patches, timesteps = timesteps,
                    start_env = c(0.3, 0.5, 0.7), temp_fluc = 0.01)

ggplot(data = e.1,
       mapping = aes(x = time, y = env1)) +
  geom_line() +
  facet_wrap(~patch) +
  theme_classic()

# simulate monocultures and mixtures in identical environmental conditions
MC <- sim_metacomm_BEF(species = species, patches = patches,
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
MC <- MC$BEF_slope

# add a identifier for the patch
MC$landscape_type <- "metacommunity"

# summarise at the landscape scale
MC_sum <- 
  MC %>%
  group_by(landscape_type, mono_mix, SR, time) %>%
  summarise(biomass = mean(biomass)) %>%
  filter(time == max(time))


# run the single patch equivalents
MC_p1 <- lapply(1:patches, function(x) {
  
  y <- sim_metacomm_BEF(species = species, patches = 1,
                        dispersal = dispersal,
                        timesteps = timesteps, 
                        start_abun = start_abun,
                        extirp_prob = extirp_prob,
                        landscape = l.1[1,], 
                        disp_mat = NA, 
                        env.df = e.1[e.1$patch == x,], 
                        env_traits.df = t.1, 
                        int_mat = si.1,
                        meas_error = 5)
  
  y <- y$BEF_slope
  
  # correct the place name
  y$patch <- x
  
  # add a identifier for the patch
  y$landscape_type <- "single_patch"
  
  return(y)
  
  })

# bind the non-interactive patches
MC_p1 <- 
  bind_rows(MC_p1) %>%
  filter(time == max(time))


MC_p1 <- 
  MC_p1 %>%
  group_by(patch) %>%
  mutate(biomass_scaled = (biomass-mean(biomass))/sd(biomass) )

split(MC_p1, MC_p1$patch) %>%
  sapply(., function(x) {
    lm.x <- lm(biomass ~ SR, data = x)
    coef(lm.x)[2]
  }) %>%
  mean(.)

MC_p1 %>%
  group_by(mono_mix, SR) %>%
  summarise(biomass = mean(biomass), .groups = "drop") %>%
  mutate(biomass_scaled = (biomass-mean(biomass))/sd(biomass) ) %>%
  lm(biomass ~ SR, data = .) %>%
  coef(.)
  
MC <- 
  MC %>%
  group_by(patch) %>%
  mutate(biomass_scaled = (biomass-mean(biomass))/sd(biomass))

split(MC, MC$patch) %>%
  sapply(., function(x) {
    lm.x <- lm(biomass ~ SR, data = x)
    coef(lm.x)[2]
  }) %>%
  mean(.)

MC_sum$biomass_scaled <- scale(MC_sum$biomass)[,1]
lm.x <- lm(biomass ~ SR, data = MC_sum)
coef(lm.x)[2]










 