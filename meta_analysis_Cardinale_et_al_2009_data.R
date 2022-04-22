
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale

# Project: Meta-analysis data from Cardinale et al. (2009)

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(ggplot2)

# make a folder to export analysis data
if(! dir.exists(here("figures"))){
  dir.create(here("figures"))
}

# link to script to draw functions from
source(here("function_plotting_theme.R"))


### download the raw data

# data come from Cardinale et al. (2009): Effects of biodiversity on the functioning of ecosystems: a summary of 164 experimental manipulations of species richness, Ecology
# link to the archive is: https://esapubs.org/archive/ecol/E090/060/#data

# load the raw data directly from the link
card.dat <- read_csv("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_1_BEF_scale_meta-analysis/BEF_summary_v2_Aug2008.csv",
                     na = c(".", "NA"),)

# check the parsing specifications
spec(card.dat)

# view the data
View(card.dat)

# subset out the relevant data points
names(card.dat)

cd.sub <- 
  card.dat %>%
  filter(FTG == "P", # primary producers
         Ygen == "SST", # standing stock measurement of ecosystem function
         Sys1 == "T", # terrestrial system
         Sys2 %in% c(3, 7, 8, 9, 10, 11, 12), # terrestrial plant system
         Des3 == "S")  # only substitutive designs

# check distribution of systems
cd.sub %>%
  group_by(Sys2) %>%
  summarise(n = n())

# based on this, let's only consider Sys2 == 7 (i.e. temperate grasslands)
cd.sub <- 
  cd.sub %>%
  filter(Sys2 == 7)


# check how many experiments this amounts to
length(unique(cd.sub$Expt))

# check range of spatial scales
range(cd.sub$SPscale, na.rm = TRUE)

# do all these data points have spatial and temporal scale data ?
cd.sub %>%
  filter(is.na(SPscale))

cd.sub %>%
  filter(is.na(Tscale))

# filter single experiment without SPscale and Tscale data
cd.sub <- 
  cd.sub %>%
  filter(!is.na(LnSPscale) | !is.na(LnTscale))

# experiments with LRR2 (where most extreme monoculture is the highest)
cd.sub <- 
  cd.sub %>%
  filter(!is.na(YEmono)) %>%
  filter(YEmono >= Y1) %>% # if the most extreme monoculture (YEmono) exceed the mean, then it is the highest
  filter(!(is.na(LRR2)))

cd.sub %>%
  select(LnSPscale, LnTscale, LRR2) %>%
  lapply(., function(x) { sum(if_else(is.na(x), 1, 0)) })


# check how many unique studies there are?
length(unique(cd.sub$Study))
length(unique(cd.sub$Expt))
nrow(cd.sub)

# how many data points come from a single experiment?
cd.sub %>%
  group_by(Study, Expt) %>%
  summarise(n = n()) %>%
  filter(n > 1)

# why is this the case? Using different time-points as individual estimates
cd.sub %>%
  filter(Study == 59, Expt == 113)

# check that this is the case for all experiments with more than one data point
cd.sub %>%
  group_by(Study, Expt) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 1) %>%
  group_by(Study, Expt) %>%
  summarise(m.t = length(unique(Dur)) )

# yes, this is the case

# for each experiment that includes multiple time points, we subset the maximum time-point (i.e. duration)
cd.sub <- 
  cd.sub %>%
  group_by(Study, Expt) %>%
  filter(Dur == max(Dur)) %>%
  ungroup()

cd.sub %>%
  group_by(Study) %>%
  summarise(n = n())

# check that there is only one experiment per data point
length(unique(cd.sub$Entry)) == length(unique(cd.sub$Expt)) 

# add raw transgressive variable to the dataset
# YSmax - response variable at the maximum species richness
# YEmono - response variable value for the most extreme monoculture
cd.sub$trans.oy <- (cd.sub$YSmax/cd.sub$YEmono)


# plot out the results: Raw

# transgressive overyielding
p1 <- 
  ggplot(data = cd.sub,
       mapping = aes(x = LnSPscale, y = (trans.oy) )) +
  geom_point() +
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  xlab("ln-spatial scale index") +
  ylab("trans. OY") +
  theme_meta()

ggsave(filename = here("figures/fig_V.pdf"), 
       plot = p1, width = 9, height = 7.5, units = "cm", dpi = 450)

cor.test(cd.sub$LnSPscale, cd.sub$trans.oy)

ggplot(data = cd.sub,
       mapping = aes(x = LnSPscale, y = LRR2)) +
  geom_point() +
  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +
  xlab("ln-spatial scale index") +
  ylab("trans. OY") +
  theme_meta()

# we also need to incorporate temporal scale as a covariate perhaps...

# next time
