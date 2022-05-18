
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale

# Project: Literature synthesis

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(ggplot2)
library(ggpubr)

# install libraries where functions are called from
if( (("gtools" %in% installed.packages()[,1]) & ("broom" %in% installed.packages()[,1]) & ("viridis" %in% installed.packages()[,1])) == FALSE) {
  print("WARNING! this script requires gtools and broom to be installed")
} 

# make a folder to export analysis data
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# link to script to draw functions from
source(here("function_plotting_theme.R"))


### download the raw data

# the raw data are archived on [Figshare](https://doi.org/10.6084/m9.figshare.12287303.v1)
# citation: Gamfeldt, Lars; Roger, Fabian; Palm, Martin; Hagan, James; Warringer, Jonas; Farewell, Anne (2020): Biodiversity and ecosystem functioning across gradients in spatial scale (literature synthesis). figshare. Dataset. https://doi.org/10.6084/m9.figshare.12287303.v1

# load the raw data
meta_dat_raw <- read_csv( url("https://ndownloader.figshare.com/files/22647539") )

# the raw data file has extra columns in the csv file structure which we remove
meta_dat_raw <- 
  meta_dat_raw %>% 
  select(-contains("...3"))

# based on our selection criteria, we removed certain data points from our original raw file

# remove 'mixture best' data points in Fridley (2003)
meta_dat_raw <- filter(meta_dat_raw, Mixture_treatment != "mixture_best")

# remove data from Fox (2002) because it is a pure resource manipulation
meta_dat_raw <- filter(meta_dat_raw, Reference != "Fox_2002")

# remove treatment manipulations that relied on disturbance
meta_dat_raw <- filter(meta_dat_raw, Env_type_1 != "disturbance")


### clean the raw data

# create a unique identifier for each experiment
meta_dat_raw <- mutate(meta_dat_raw, Experiment_ID = paste(Reference, Experiment_number, sep = "_"))

# translate ecosystem function values to positive if low value indicates high ecosystem function
meta_dat_raw <- 
  meta_dat_raw %>% 
  group_by(Experiment_ID) %>% 
  mutate(ef_min = min(Ecosystem_function_mean)) %>% 
  ungroup() %>%
  mutate(ef_min = if_else(ef_min < 0, (-ef_min), 0)) %>% 
  mutate(Ecosystem_function_mean = (Ecosystem_function_mean + ef_min)) %>%
  select(-ef_min)


### prepare a copy of the data for analysis

# create a copy for the analysis with reordered columns
meta_an <- 
  meta_dat_raw %>%
  select(Experiment_ID, names(meta_dat_raw))


### analysis

# for each Experiment_ID, calculate log(mixture/mean monoculture) for each combination of environments

# convert to list by experiment ID
meta_l <- split(meta_an, meta_an$Experiment_ID)

l.out <- 
  lapply(meta_l, function(r) {
    
    d.l.in <- r
    
    # get all combinations of environments
    x <- unique(d.l.in$Environment)
    c.n <- vector("list", length = length(x))
    for (j in x) {
      c.n[[j]] <- gtools::combinations(n = length(x), r = j, v = x, repeats.allowed = FALSE)
    }
    
    
    r.l <- lapply(c.n, function(p) {
      
      nbe.i <- vector("list", length = nrow(p))
      
      for (i in 1:nrow(p)) {
        
        # get rows that are relevant
        n.r <- which(d.l.in$Environment %in% p[i,])
        
        # get relevant data.frame
        d.f <- d.l.in[n.r, ]
        
        # summarise data using dplyr
        d.f.s <- 
          d.f %>% 
          group_by(Mixture_treatment, Mixture_ID) %>%
          summarise(m.func = mean(Ecosystem_function_mean, na.rm = TRUE), .groups = "drop")
        
        # calculate ln response ratio of mixture and transgressive overyielding
        
        # get mixture
        m.x <- d.f.s[d.f.s$Mixture_treatment == "mixture", ]$m.func
        
        # get mean monoculture
        m.m <- mean(d.f.s[d.f.s$Mixture_treatment == "monoculture", ]$m.func, na.rm = TRUE)
        
        # get best monoculture
        m.b <- max(d.f.s[d.f.s$Mixture_treatment == "monoculture", ]$m.func, na.rm = TRUE)
        
        # get ln-response ratio
        ln.rr <- log((m.x/m.m))
        
        # get transgressive overyielding
        t.oy <- (m.x/m.b)
        
        # combine these into a vector
        nbe.i[[i]] <- c(ln.rr, t.oy)
        
      }
      
      return(nbe.i)
      
    })
    
    # bring this into a more manageable format
    
    # unlist the environmental combinations
    l.s <- unlist(lapply(c.n, function(x) {apply(x, 1, paste, collapse="")}))
    
    # unlist the lrr and toy
    r.s <- unlist(r.l)
    
    # pull this into a data.frame
    d.f.o <- data.frame(environment.id = rep(l.s, each = 2),
                        response = rep(c("l.rr", "t.oy"), times = length(l.s)),
                        value = r.s)
    
    # convert l.rr and t.oy into separate columns using pivot_wider()
    d.f.o <- 
      pivot_wider(d.f.o,
                  names_from = "response",
                  values_from = "value")
    
    # add a column for the number of environments
    d.f.o <- 
      d.f.o %>%
      mutate(environments = nchar( as.character(environment.id) ))
    
    # reorder the columns
    d.f.o <- 
      d.f.o %>%
      select(environment.id, environments, l.rr, t.oy)
    
    # return this data.frame
    return(d.f.o)
    
  })


### bind this into a data.frame by Experiment_ID
l.df <- bind_rows(l.out, .id = "Experiment_ID")


### calculate the slope between scale and (i) l.rr, and (ii) t.oy for each Experiment_ID

# get slope between scale and l.rr and t.oy for each Experiment_ID
l.est <- 
  lapply(l.out, function(x) {
    
    lm.x <- lm(l.rr ~ environments, data = x)
    lm.y <- lm(t.oy ~ environments, data = x)
    
    data.frame(l.rr.est = as.numeric(lm.x$coefficients[2]),
               t.oy.est = as.numeric(lm.y$coefficients[2]))
    
  })

# bind these slope estimates into a data.frame
l.est <- bind_rows(l.est, .id = "Experiment_ID")


### calculate exploratory covariates

# create a data.frame from which to calculate the covariates
cov.df.prep <- 
  meta_an %>%
  select(Experiment_ID, Environment, Mixture_treatment, Mixture_ID, Ecosystem_function_mean)

# calculate average and coefficient of variation in average overyielding across environments
m.oy <- 
  cov.df.prep %>%
  group_by(Experiment_ID, Environment, Mixture_treatment) %>%
  summarise(function.m = mean(Ecosystem_function_mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = "Mixture_treatment", values_from = "function.m") %>%
  mutate(overyielding = mixture/monoculture) %>%
  group_by(Experiment_ID) %>%
  summarise(overyielding.m = mean(overyielding, na.rm = TRUE))

# calculate the species specialisation index
ss.i <- 
  cov.df.prep %>%
  filter(Mixture_treatment == "monoculture") %>%
  group_by(Experiment_ID, Environment) %>%
  filter(Ecosystem_function_mean == max(Ecosystem_function_mean)) %>%
  ungroup() %>%
  group_by(Experiment_ID) %>%
  summarise(n.spp = length(unique(Mixture_ID)),
            n.env = length(unique(Environment)), .groups = "drop") %>%
  mutate(n.spp.prop = (n.spp/n.env) ) %>% 
  mutate(species.specialisation = (n.spp.prop - (1/n.env))/(1-(1/n.env)) ) %>%
  select(Experiment_ID, species.specialisation)

# join these data.frames   
cov.df <- 
  full_join(m.oy, ss.i, by = "Experiment_ID")


### join the covariates to the slopes between scale and l.rr and t.oy
est.cov <- 
  full_join(cov.df, l.est, by = "Experiment_ID")


### plots and statistics for manuscript

x1 <- 
  ggplot(data = l.df,
       mapping = aes(x = environments, y = t.oy, colour = Experiment_ID)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, size = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.2) +
  ylab("Trangressive overyielding") +
  xlab("Scale") +
  scale_colour_viridis_d(option = "C") +
  theme(legend.position = "none") +
  theme_meta()

x2 <- 
  ggplot(data = l.est,
       mapping = aes(x = t.oy.est)) +
  ylab("") +
  xlab("Transgressive overyielding ~ scale (est.)") +
  geom_histogram(alpha = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  # scale_x_continuous(limits = c(-0.5, 0.6)) +
  theme_meta()

# plot raw correlation between average overyielding and species specialisation
x3 <- 
  ggplot(data = est.cov,
       mapping = aes(x = species.specialisation, y = t.oy.est, colour = overyielding.m)) +
  geom_jitter(width = 0.01) +
  ylab("Trangressive overyielding ~ scale (est.)") +
  xlab("Species specialisation index") +
  viridis::scale_colour_viridis(option = "C", end = 0.9) +
  guides(color = guide_colourbar(title.position = "top", 
                                 title.vjust = 1,
                                 frame.colour = "black", 
                                 ticks.colour = NA,
                                 barwidth = 5,
                                 barheight = 0.3)) +
  labs(col="Average overyielding") +
  theme_meta() +
  theme(legend.position = c(0.6, 0.75),
        legend.direction="horizontal",
        legend.justification=c(1, 0), 
        legend.key.width=unit(1, "lines"), 
        legend.key.height=unit(1, "lines"),
        # plot.margin = unit(c(3, 1, 0.5, 0.5), "lines"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

# plot raw correlation between average overyielding and species specialisation
x4 <- 
  ggplot(data = est.cov,
       mapping = aes(x = overyielding.m, y = t.oy.est, colour = species.specialisation)) +
  geom_point() +
  ylab("Trangressive overyielding ~ scale (est.)") +
  xlab("Average overyielding") +
  viridis::scale_colour_viridis(option = "C", end = 0.9) +
  guides(color = guide_colourbar(title.position = "top", 
                                 title.vjust = 1,
                                 frame.colour = "black", 
                                 ticks.colour = NA,
                                 barwidth = 5,
                                 barheight = 0.3)) +
  labs(col="Species specialisation") +
  theme_meta() +
  theme(legend.position = c(0.6, 0.75),
        legend.direction="horizontal",
        legend.justification=c(1, 0), 
        legend.key.width=unit(1, "lines"), 
        legend.key.height=unit(1, "lines"),
        # plot.margin = unit(c(3, 1, 0.5, 0.5), "lines"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))

fY <- ggpubr::ggarrange(x1, x2, x3,x4, ncol = 2, nrow = 2,
                  labels = c("A", "B", "C", "D"),
                  font.label = list(size = 10, color = "black", face = "bold", family = NULL))

ggsave(filename = here("figures/Fig_4.pdf"), 
       plot = fY, width = 18, height = 16, units = "cm", dpi = 450)

# perform wilcoxon tests on scale-slopes

# defines slopes to test
s.t <- c("t.oy.est")

# define alternative hypotheses
a.h <- c("two.sided")

wc.test <- vector("list")
for(i in 1:length(s.t) ) {
  
  w <- wilcox.test(est.cov[[ s.t[i] ]], alternative = a.h[i], 
                   mu = 0, paired = FALSE, exact = FALSE, conf.int = FALSE)
  wc.test[[i]] <- 
    data.frame(response = s.t[i], 
               w.stat = w$statistic,
               p.value = w$p.value, 
               row.names = NULL)
}
wc.test <- do.call(rbind, wc.test)

# correct the p.values for multiple testing using bonferroni
wc.test$p.value.bf <- p.adjust(wc.test$p.value, method = "bonferroni")

# is the trangressive overyielding result robust to outliers?
est.cov %>%
  filter(t.oy.est != max(t.oy.est)) %>%
  pull(t.oy.est) %>%
  wilcox.test(., alternative = "two.sided", 
              mu = 0, paired = FALSE, exact = FALSE, conf.int = FALSE)


### plot fig. S5

# add information on realm and function-type
si.dat <- 
  full_join(l.df,
            meta_an %>%
              select(Experiment_ID, Realm, `Ecosystem function` = Ecosystem_function_type) %>%
              distinct(),
            by = "Experiment_ID")

# rename the experiment_ID column entries
si.dat$Experiment_ID <- factor(si.dat$Experiment_ID)
levels(si.dat$Experiment_ID) <- LETTERS[1:length(si.dat$Experiment_ID)]

f.S5 <- 
  ggplot(data = si.dat, 
         mapping = aes(x = environments, y = t.oy,
                       colour = Realm, shape = `Ecosystem function`,
                       group = NULL)) +
  geom_jitter(width = 0.01, size = 1.2, alpha = 0.7) +
  ylab("Transgressive overyielding") +
  xlab("Scale") +
  scale_colour_viridis_d(option = "C", end = 0.9) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Experiment_ID, scales = "free", ncol = 5) +
  labs(shape = "Ecosystem function", colour = "Realm") +
  scale_x_continuous(breaks = c(1:8)) +
  theme_meta() +
  theme(axis.text = element_text(size = 6),
        legend.position = c(0.8, 0.045),
        legend.box = "horizontal",
        legend.key = element_rect(fill = "white"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 8))
f.S5

# save the output
ggsave(filename = here("figures/Fig_S5.pdf"), 
       plot = f.S5, width = 15.24, height = 20.32, units = "cm", dpi = 450)

### END
