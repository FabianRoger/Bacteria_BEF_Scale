
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale

# Project: Literature synthesis (meta-analysis)

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(here)

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
  select(-contains("X3"))

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
      mutate(environments = nchar(environment.id))
    
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

# calculate coefficient of variation among environments in mean function across monocultures and mixtures
f.cv <- 
  cov.df.prep %>%
  group_by(Experiment_ID, Environment) %>%
  summarise(function.m = mean(Ecosystem_function_mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(Experiment_ID) %>%
  summarise(function.cv = (sd(function.m)/mean(function.m))*100, .groups = "drop")

# calculate average and coefficient of variation in average overyielding across environments
m.cv.oy <- 
  cov.df.prep %>%
  group_by(Experiment_ID, Environment, Mixture_treatment) %>%
  summarise(function.m = mean(Ecosystem_function_mean, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = "Mixture_treatment", values_from = "function.m") %>%
  mutate(overyielding = mixture/monoculture) %>%
  group_by(Experiment_ID) %>%
  summarise(overyielding.m = mean(overyielding, na.rm = TRUE),
            overyielding.cv = (sd(overyielding)/mean(overyielding))*100, .groups = "drop")

# join these data.frames
cov.df <- full_join(f.cv, m.cv.oy, by = "Experiment_ID")

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
  full_join(cov.df, ss.i, by = "Experiment_ID")


### join the covariates to the slopes between scale and l.rr and t.oy
est.cov <- 
  full_join(cov.df, l.est, by = "Experiment_ID")


### plots and statistics for the manuscript

# perform wilcoxon tests on scale-slopes

# defines slopes to test
s.t <- c("l.rr.est", "t.oy.est")

# define alternative hypotheses
a.h <- c("two.sided", "two.sided")

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


# load the ggplot and ggpubr libraries for plotting
library(ggplot2)
library(ggpubr)
library(viridis)


# Ecology (journal) figure guidelines:
# (1) portrait layout (maximum 6 inches (15.24 cm) wide x 8 inches (20.32 cm) high)
# (2) landscape layout (maximum 8.75 inches (22.225 cm) wide x 5.25 inches (13.335 cm) high)


### fig. 3

# choose example Experiment_ID
e.g.s <- "Blake_and_Duffy_2010_one"
e.g.s.l <- "Blake and Duffy 2010"

# get a middle viridis colour
v.c <- viridis(n = 1, alpha = 1, begin = 0.5, end = 0.5, option = "C")
v.c.1 <- viridis(n = 1, alpha = 1, begin = 0.7, end = 0.7, option = "C")

# fig. 3a
f.3a <- 
  ggplot(data = filter(est.cov, Experiment_ID != e.g.s),
       mapping = aes(x = 0.9, y = l.rr.est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  geom_point(size = 2, position = position_jitter(w = 0.02, h = 0), 
             shape = 16, alpha = 0.5) +
  scale_x_continuous(limits = c(0.85, 1.05)) +
  geom_point(mapping = aes(x = 0.91, y = est.cov[est.cov$Experiment_ID == e.g.s, ]$l.rr.est ), 
             colour = v.c, size = 2, shape = 16) +
  # geom_point(mapping = aes(x = 0.9, y = median(est.cov$l.rr.est)),
             # shape = 21, size = 3.5, fill = v.c.1, colour = "black") +
  geom_hline(yintercept = median(est.cov$l.rr.est), size = 0.5, colour = v.c.1) +
  annotate(geom = "text", x = 1.025, y = 0.25, 
           label = paste("W", wc.test[wc.test$response == "l.rr.est", ]$w.stat, sep = " = "),
           size = 2.5) +
  annotate(geom = "text", x = 1.025, y = 0.2, 
           label = paste("P", round(wc.test[wc.test$response == "l.rr.est", ]$p.value.bf, 3), sep = " = "),
           size = 2.5) +
  ylab("log(mix/mono) ~ scale (est.)") +
  xlab("") +
  ggtitle(label = "") +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_blank())

# fig. 3b
f.3b <- 
  l.df %>%
  filter(Experiment_ID == e.g.s ) %>%
  ggplot(data = .,
         mapping = aes(x = environments, y = l.rr)) +
  geom_smooth(method = "lm", se = FALSE, colour = v.c, size = 0.75) +
  geom_jitter(width = 0.01, size = 2, 
              colour = v.c, shape = 16) +
  xlab("scale") +
  ylab("log(mix/mono)") +
  ggtitle(label = e.g.s.l) +
  theme_meta()

# fig. 3c
f.3c <- 
  ggplot(data = filter(est.cov, Experiment_ID != e.g.s),
       mapping = aes(x = 0.9, y = t.oy.est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25) +
  geom_point(size = 2, position = position_jitter(w = 0.02, h = 0), 
              shape = 16, alpha = 0.5) +
  scale_x_continuous(limits = c(0.85, 1.05)) +
  scale_y_continuous(limits = c(-0.1, 0.5)) +
  geom_point(mapping = aes(x = 0.91, y = est.cov[est.cov$Experiment_ID == e.g.s, ]$t.oy.est ), 
             colour = v.c, size = 2, shape = 16) +
  # geom_point(mapping = aes(x = 0.9, y = median(est.cov$t.oy.est)),
             # shape = 21, size = 3.5, fill = v.c.1, colour = "black") +
  geom_hline(yintercept = median(est.cov$t.oy.est), size = 0.5, colour = v.c.1) +
  annotate(geom = "text", x = 1.025, y = 0.45, 
           label = paste("W", wc.test[wc.test$response == "t.oy.est", ]$w.stat, sep = " = "),
           size = 2.5) +
  annotate(geom = "text", x = 1.025, y = 0.4, 
           label = paste("P", round(wc.test[wc.test$response == "t.oy.est", ]$p.value.bf, 3), sep = " = "),
           size = 2.5) +
  ylab("trans. OY ~ scale (est.)") +
  xlab("") +
  ggtitle(label = "") +
  theme_meta() +
  theme(axis.text.x = element_text(colour = "white"),
        axis.ticks.x = element_blank())

# fig. 3d
f.3d <- 
  l.df %>%
  filter(Experiment_ID == e.g.s) %>%
  ggplot(data = .,
         mapping = aes(x = environments, y = t.oy)) +
  geom_smooth(method = "lm", se = FALSE, colour = v.c, size = 0.75) +
  geom_jitter(width = 0.01, size = 2, 
              colour = v.c, shape = 16) +
  xlab("scale") +
  ylab("trans. OY") +
  ggtitle(label = e.g.s.l) +
  theme_meta()

f.3 <- 
  ggarrange(f.3a, f.3b, f.3c, f.3d, ncol = 2, nrow = 2,
            widths = c(1.3, 1) )

ggsave(filename = here("figures/fig_3.jpg"), 
       plot = f.3, width = 11, height = 11, units = "cm")


# is the trangressive overyielding result robust to outliers?
est.cov %>%
  filter(t.oy.est != max(t.oy.est)) %>%
  pull(t.oy.est) %>%
  wilcox.test(., alternative = "two.sided", 
              mu = 0, paired = FALSE, exact = FALSE, conf.int = FALSE)


### how is is t.oy.est related to covariates?

# set up a function to run different models that can then be compared
lm.scale <- function(data, slope, e.vars, outliers = NA) {
  
  # set an output list for the model coefficients
  est.lm <- vector("list", length(e.vars))
  names(est.lm) <- seq(1:length(e.vars))
  
  # set an output list for the model fit statistics
  fit.lm <- vector("list", length(e.vars))
  names(fit.lm) <- seq_along(1:length(e.vars))
  
  for (i in 1:length(e.vars) ) {
    
    # subset out outliers
    df <- 
      data %>%
      filter(!Experiment_ID %in% outliers)
    
    # fit model using chosen predictors
    lm.e.vars <- lm(reformulate(e.vars[[i]], slope), data = df)
    
    # write coefficients to the est.lm list
    est.lm[[i]] <- broom::tidy(lm.e.vars)
    
    # write fit statistics to the fit.lm list
    fit.lm[[i]] <- broom::glance(lm.e.vars)
  }
  
  # convert lists to data.frames and join
  full_join(bind_rows(est.lm, .id = "model"), 
            bind_rows(fit.lm, .id = "model"),
            by = "model")
}

# set up the explanatory variables for the different models
exp.vars <- list(c("species.specialisation", "overyielding.m", "function.cv", "overyielding.cv"),
                 c("species.specialisation*overyielding.m"), 
                 c("species.specialisation*function.cv"),
                 c("species.specialisation*overyielding.cv"),
                 c("species.specialisation"),
                 c("overyielding.cv"),
                 c("function.cv"),
                 c("overyielding.m"))

# run this set of models with all the data
t.oy.exp <- lm.scale(data = est.cov, slope = "t.oy.est", e.vars = exp.vars, outliers = NA)

# clean this output
table.s3 <- 
  t.oy.exp %>% 
  select(model, term, r.squared, AIC) %>%
  group_by(model, r.squared, AIC) %>%
  summarise(terms = paste(term, collapse = "+")) %>% 
  ungroup() %>%
  arrange(AIC) %>% 
  select(terms, r.squared, AIC) %>%
  mutate(delta_AIC =  AIC - (min(AIC))) %>%
  mutate(AIC_wt_start = exp(1)^(-0.5*delta_AIC)) %>%
  mutate(AIC_wt = AIC_wt_start/sum(AIC_wt_start)) %>%
  select(-AIC_wt_start)

# write_csv(table.s3, here("figures/Table_S3.csv") )


# re-run this analysis without the large outlier
t.oy.exp.outlier <- lm.scale(data = est.cov, slope = "t.oy.est", e.vars = exp.vars, outliers = c("Dzialowski_Smith_2008_one"))

# clean this output
table.s4 <- 
  t.oy.exp.outlier %>% 
  select(model, term, r.squared, AIC) %>%
  group_by(model, r.squared, AIC) %>%
  summarise(terms = paste(term, collapse = "+")) %>% 
  ungroup() %>%
  arrange(AIC) %>% 
  select(terms, r.squared, AIC) %>%
  mutate(delta_AIC =  AIC - (min(AIC))) %>%
  mutate(AIC_wt_start = exp(1)^(-0.5*delta_AIC)) %>%
  mutate(AIC_wt = AIC_wt_start/sum(AIC_wt_start)) %>%
  select(-AIC_wt_start)

# write_csv(table.s4, here("figures/Table_S4.csv") )


### plot fig.4

# fit the linear model with the lowest AIC (highest AIC weight)
t.oy.exp %>%
  filter(AIC == min(AIC)) %>%
  pull(term)

lm.b <- lm(t.oy.est ~ species.specialisation*overyielding.m, data = est.cov)

# check model predictions
plot(lm.b)

# get predicted values from this model

# species.specialisation gradient
range(est.cov$species.specialisation)
ssi.g <- 
  seq(from = min(est.cov$species.specialisation), 
      to = max(est.cov$species.specialisation), 
      by = 0.1)

# set three levels of overyielding.m
range(est.cov$overyielding.m)
oy.m.g <- c(0.75, 1, 2)

# put these values into a data.frame
g.pred <- expand.grid(species.specialisation = ssi.g, overyielding.m = oy.m.g)

# get the predicted values
g.pred$t.oy.est <- predict(object = lm.b, newdata = g.pred)


# plot fig. 4
f.4 <- 
  est.cov %>%
  rename("average overyielding" = overyielding.m) %>%
  ggplot(data = .,
         mapping = aes(x = species.specialisation, y = t.oy.est, colour = `average overyielding`)) +
  geom_line(data = rename(g.pred, "average overyielding" = overyielding.m), 
            aes(group = (`average overyielding`)),
            size = 0.75) +
  geom_jitter(width = 0.01, size = 3) +
  ylab("trans. OY ~ scale (est.)") +
  xlab("species specialisation index") +
  viridis::scale_colour_viridis() +
  annotate(geom = "text", x = 0.1, y = 0.35, size = 3.5,
           label = paste("r^2 == ", round(summary(lm.b)$r.squared, 2)), parse = TRUE) +
  guides(color = guide_colourbar(title.position = "top", 
                                 title.vjust = 1,
                                 frame.colour = "black", 
                                 ticks.colour = NA,
                                 barwidth = 10,
                                 barheight = 0.5)) +
  theme_meta() +
  theme(legend.position = c(0.6, 0.75),
        legend.direction="horizontal",
        legend.justification=c(1, 0), 
        legend.key.width=unit(1, "lines"), 
        legend.key.height=unit(1, "lines"),
        # plot.margin = unit(c(3, 1, 0.5, 0.5), "lines"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

ggsave(filename = here("figures/fig_4.jpg"), 
       plot = f.4, width = 14, height = 12, units = "cm")





