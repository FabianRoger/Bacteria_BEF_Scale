
# Title: Gamfeldt et al. Biodiversity and ecosystem functioning across gradients in spatial scale

# Project: Meta-analysis data assessment

# load relevant libraries
install_if <- function(x) {
  
  if(x %in% rownames(installed.packages()) == FALSE) {
    message(paste(x, "is equired but not installed. Installing now"))
    Sys.sleep(1)
    install.packages(x)
    library(x)
  } else{ 
    library(x, character.only=TRUE)}
}

install_if("dplyr")
install_if("readr")


### download data

# the raw data are archived on [Figshare](https://doi.org/10.6084/m9.figshare.12287303.v1)
# citation: Gamfeldt, Lars; Roger, Fabian; Palm, Martin; Hagan, James; Warringer, Jonas; Farewell, Anne (2020): Biodiversity and ecosystem functioning across gradients in spatial scale (literature synthesis). figshare. Dataset. https://doi.org/10.6084/m9.figshare.12287303.v1

# load the raw data
data_ass <- read_csv( url("https://ndownloader.figshare.com/files/22647536") )

# check the data
View(data_ass)

# check the variable names
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

### END

