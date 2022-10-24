# Script initialised 2018-10-30
# Lars Gamfeldt

# Sorting and exploring bacteria straindata for BEF scale experiment. 
# Data are for screening strains for choosing which strains to use in experiment. 
# Data file is for BSL2_6_7_ATCC, absolute data. 
# 4 antibiotics: cam, kan, tet, trim. 2 plates for each.
# 1 control, 2 plates. 
# For each, the file includes yield (e.g. yield_antibiotic1) and generation time (e.g. gentime_antibiotic1). 

# Libraries

install_if <- function(x) {
  
  if(x %in% rownames(installed.packages()) == FALSE) {
    message(paste(x, "is required but not installed. Installing now"))
    Sys.sleep(1)
    install.packages(x)
    library(x)
  } else{ 
    library(x, character.only=TRUE)}
}


install_if("dplyr")
install_if("tidyr")
install_if("readr")
install_if("ggplot2")
install_if("tibble")
install_if("vegan")



# Read in data
rawdat<- read_delim("Data/BSL2_6_7_all_ABs_absolute.csv", delim = ";")
compiled_cam_kan<- read_delim("Data/181019_BSL2-6-7_Cam-1_Cam-1_Kan-8_Kan-8_M9_Cersei_compiled.csv", delim = ",")
compiled_tet_trim<- read_delim("Data/181019_BSL2-6-7_Tet-0.3_Tet-0.3_Trim-0.25_Trim-0.25_M9_Jamie_compiled.csv", delim = ",")
compiled_con<- read_delim("Data/181019_BSL2-6-7_Con_Con_Non_Non_M9_Tyrion_compiled.csv", delim = ",")

#extract 
dat<- bind_rows(compiled_cam_kan, compiled_tet_trim)
dat$Plate<- as.character(dat$Plate)
dat$Position<- as.character(gsub("^._","",dat$Position))


#check for outliers
dat %>% 
  group_by(Position, Strain, Condition, Plate) %>% 
  summarise(max_value = max(value, na.rm = TRUE),
            max_quant = quantile(value, 0.8, na.rm = TRUE)) %>% 
  ggplot(aes (y = max_value, x = max_quant))+
  geom_point()




# Retrieve max value for each strain
dat<- dat %>% 
  group_by(Position, Strain, Condition, Plate) %>% 
  summarise(max_value = max(value, na.rm = TRUE))


# exclude strains
dat <- dat %>%
  ungroup() %>% 
  filter(!grepl("ATCC", Strain), !grepl("EMPTY", Strain)) 

# visualize difference between plates
dat %>% 
  ggplot(aes(x = Plate, y = max_value, group = Strain))+
  geom_point()+
  geom_line(alpha = 0.3)+
  scale_y_log10()+
  facet_wrap(~Condition, scales = "free")+
  theme_bw()

# check differences between plates
dat_filt <- 
dat %>% 
  mutate(Plate = case_when(Plate == "0" ~ "0",
                           Plate == "1" ~ "1",
                           Plate == "2" ~ "0",
                           Plate == "3" ~ "1",)) %>%
  spread(Plate, max_value) %>% 
  mutate(diff = `1` / `0`) %>% 
  filter(`0` != -Inf & `1` != -Inf) %>% 
  mutate(diff = ifelse(diff<1, 1/diff, diff)) 

dat_filt %>% 
  filter(diff < 10) %>% 
  group_by(Condition) %>% 
  mutate(rank_diff = rank(diff)) %>% 
  ggplot(aes(x = diff, y = rank_diff))+
  facet_wrap(~Condition)+
  geom_point()

#filter out all strains with differences larger than a factor of 2.5 between replicate plates
#calculate mean value between plates

dat_filt <- 
dat_filt %>% 
  filter(diff <= 2.5) %>% 
  mutate(max_value = (`0`+`1`)/2) %>% 
  select(Strain, Condition, max_value)
  
#keep only strains that perform in the middle 50% on controle plates
intermediate_strains <- 
compiled_con %>% 
  group_by(Plate, Strain) %>% 
  filter(!grepl("EMPTY", Strain),
         !grepl("ATCC", Strain)) %>% 
  summarise(max_value = max(value)) %>% 
  group_by(Strain) %>%
  filter(!any(is.na(max_value))) %>% 
  group_by(Strain) %>% 
  summarise(max_value = mean(max_value)) %>% 
  filter(max_value > quantile(max_value, 0.25) &
           max_value > quantile(max_value, 0.75))

dat_filt <- 
dat_filt %>% 
  filter(Strain %in% intermediate_strains$Strain)

# # Convert columns from factor to numeric
# dat <- 
#   rawdat %>% 
#   mutate_at(vars(contains("yield"), contains("gentime")), 
#             function(x) as.numeric(as.character(x))) 
# 
# # Drop all rows that contain ATCC
# dat<- dat[- grep("ATCC", dat$strain),]
# dat<- dat[- grep("EMPTY", dat$strain),]
# 
# # Add variables with mean yield and generation time for each antibiotic and control
# # For yield
# 
# dat_long <- 
# dat %>% 
#   gather(variable, value, -strain) %>% 
#   mutate(var = gsub("(\\w+)_.+", "\\1", variable)) %>% 
#   mutate(antibio = gsub("\\w+_(\\w+)\\d", "\\1", variable)) %>% 
#   mutate(rep = gsub("\\w+_\\w+(\\d)","\\1", variable)) %>% 
#   select(-variable)
# 
# #remove all strains with any NA values
# dat_long <- 
#   dat_long %>% 
#   group_by(strain) %>% 
#   filter(!any(is.na(value)))
# 
# #calculate mean of replicates
# dat_long<- dat_long %>% 
#   group_by(strain, var, antibio) %>% 
#   summarise(value = mean(value))
# 
# #plot the distribution of values
# dat_long %>% 
#   filter(var == "gentime") %>% 
#   ggplot(aes(x = value))+
#   geom_histogram()+
#   facet_wrap(~antibio)+
#   theme_bw()+
#   labs(title = "gentime")
# 
# dat_long %>% 
#   filter(var == "yield") %>% 
#   ggplot(aes(x = value))+
#   geom_histogram()+
#   facet_wrap(~antibio)+
#   theme_bw()+
#   scale_x_log10()+
#   labs(title = "yield")
# 
# #exclude controle
# dat_long <- dat_long %>% 
#   filter(antibio != "con")



# calculate cv of yield
max100cv <- 
  dat_filt %>% 
  group_by(Strain) %>% 
  filter(n()==4) %>% #314 strains fulfill this requirement
  filter(!any(max_value < 2e8)) %>% #84 Strains fulfill this requirement (maybe too strict?)
  summarise(
    cv = sd(max_value)/mean(max_value),
    value = mean(max_value)) %>% 
  ungroup()


dat_cv <- 
dat_filt %>% 
  #filter(var == "yield") %>% 
  filter(Strain %in% max100cv$Strain) %>% 
  ungroup() %>% 
 # select(Strain, Condition, max_value) %>% 
  spread(Condition, max_value) 

dat_cv_mat <- as.matrix(dat_cv[, -1]) 

dimnames(dat_cv_mat) <-  list(
  as.character(dat_cv$Strain),
  colnames(dat_cv[,-1])) 

#standardize by row total
dat_cv_mat <- decostand(dat_cv_mat, margin = 1, method = "total")

yield_dist <- vegdist(dat_cv_mat, method = "eu") %>% 
  as.matrix()

i <- 1
yield_dist_set <- yield_dist
max_names <- c()

while(i < 6){
  i <- i+1
  temp_names <- rownames(which(yield_dist_set == max(yield_dist_set), arr.ind = T))
  max_names <- c(max_names, temp_names)
  yield_dist_set <- yield_dist[!colnames(yield_dist) %in% max_names, max_names]
}

NMDS <- metaMDS(dat_cv_mat, dist = "eu")
scrs <- as.data.frame(scores(NMDS, display = "sites"))

vf <- envfit(NMDS, dat_cv_mat, perm = 999)
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))

ggplot(scrs) +
  geom_point(mapping = aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() + ## need aspect ratio of 1!
  geom_segment(data = spp.scrs,
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")), colour = "grey") +
  geom_text(data = spp.scrs, aes(x = NMDS1, y = NMDS2, label = Species),
            size = 3)

# 
# 
# top_50 <- 
#   apply(yield_dist, 2, sum) %>% 
#   sort(decreasing = T)
# 
# top_50 <- top_50[1:50]
# 
# spec_comb <- combn(names(top_50), 6)
#   
# max_dist <- apply(spec_comb, 2, function(x) sum(yield_dist[x,x]))
# 
# max_dist_df <- 
# tibble(comb = apply(spec_comb, 2, list),
#        dist = max_dist) %>% 
#   arrange(desc(dist))
# 
# top1 <- unlist(max_dist_df[1,1])

dat_filt %>% 
  filter(Strain %in% max_names) %>% 
  ggplot(aes(x = Condition, y = max_value, fill = Condition))+
  geom_bar(stat = "identity")+
  facet_wrap(~Strain)+
  theme_bw()+
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "Set1")
  

# ################################
# 
# dat<- dat %>% 
#   mutate(yield_kan = (yield_kan1 + yield_kan2)/2)
# 
# dat<- dat %>% 
#   mutate(yield_tet = (yield_tet1 + yield_tet2)/2)
# 
# dat<- dat %>% 
#   mutate(yield_trim = (yield_trim1 + yield_trim2)/2)
# 
# dat<- dat %>% 
#   mutate(yield_con = (yield_con1 + yield_con2)/2)
#        
# # For generation time
# dat<- dat %>% 
#   mutate(gentime_cam = (gentime_cam1 + gentime_cam2)/2)
# dat<- dat %>% 
#   mutate(gentime_kan = (gentime_kan1 + gentime_kan2)/2)
# dat<- dat %>% 
#   mutate(gentime_tet = (gentime_tet1 + gentime_tet2)/2)
# dat<- dat %>% 
#   mutate(gentime_trim = (gentime_trim1 + gentime_trim2)/2)
# dat<- dat %>% 
#   mutate(gentime_con = (gentime_con1 + gentime_con2)/2)
# 
# # Inspect dat
# str(dat)
# summary(dat)
# 
# # Inspect histogram of yield
# hist(log(dat$yield_con), breaks = 50)
# 
# # First step is to subset data by yield in control (only glucose 0.1%)
# # Subsetting by > median and < 90% quantile. 
# median(dat$yield_con, na.rm = TRUE) 
# median(log(dat$yield_con), na.rm = TRUE) # median of log yield
# mean(log(dat$yield_con), na.rm = TRUE) # mean of log yield
# 
# dat1<- filter(dat, yield_con > median(yield_con, na.rm = TRUE) & yield_con < quantile(yield_con, 0.9, na.rm = TRUE))
# dat1b<- filter(dat, log(yield_con) > quantile(log(yield_con), 0.2, na.rm = TRUE) & log(yield_con) < quantile(log(yield_con), 0.8, na.rm = TRUE))
# 
# dim(dat1b)
# hist(dat1b$yield_con, breaks = 50)
# summary(dat1b)
# 
# # Select subset of colums
# dat2<- select(dat1b, c(strain, yield_cam:yield_con, gentime_cam:gentime_con))
# 
# # Remove rows with NA
# dat3<- dat2 %>% 
#   na.omit()
# 
# # Inspect dat3
# str(dat3)
# summary(dat3)
# 
# # Inspect histogram of yield
# hist(dat3$yield_con, breaks = 50)
# 
# # Sorting (arranging) the data based on top 50% in each antibiotic
# dat_cam<- dat3 %>% 
#   arrange(desc(yield_cam))
# dat_kan<- dat3 %>% 
#   arrange(desc(yield_kan))
# dat_tet<- dat3 %>% 
#   arrange(desc(yield_tet))
# dat_trim<- dat3 %>% 
#   arrange(desc(yield_trim))
# 
# dat_cam<- add_column(dat_cam, sorted = "cam")
# dat_kan<- add_column(dat_kan, sorted = "kan")
# dat_tet<- add_column(dat_tet, sorted = "tet")
# dat_trim<- add_column(dat_trim, sorted = "trim")
# 
# View(dat_cam)
# View(dat_kan)
# View(dat_tet)
# View(dat_trim)
# 
# # Filter strains for each antibiotic that is in top 70%, and then above bottom 20% and less than 60% for the other three antibiotics
# dat_cam_best<- filter(dat_cam, yield_cam > quantile(yield_cam, 0.7) & 
#                     yield_kan < quantile(yield_kan, 0.6) & yield_kan > quantile(yield_kan, 0.2) &
#                     yield_tet < quantile(yield_tet, 0.6) & yield_tet > quantile(yield_tet, 0.2) &
#                     yield_trim < quantile(yield_trim, 0.6) & yield_trim > quantile(yield_trim, 0.2))
# dim(dat_cam_best)
# 
# dat_kan_best<- filter(dat_kan, yield_kan > quantile(yield_kan, 0.7) & 
#                         yield_cam < quantile(yield_cam, 0.6) & yield_cam > quantile(yield_cam, 0.2) &
#                         yield_tet < quantile(yield_tet, 0.6) & yield_tet > quantile(yield_tet, 0.2) &
#                         yield_trim < quantile(yield_trim, 0.6) & yield_trim > quantile(yield_trim, 0.2))
# dim(dat_kan_best)
# 
# dat_tet_best<- filter(dat_tet, yield_tet > quantile(yield_tet, 0.7) & 
#                         yield_kan < quantile(yield_kan, 0.6) & yield_kan > quantile(yield_kan, 0.2) &
#                         yield_cam < quantile(yield_cam, 0.6) & yield_cam > quantile(yield_cam, 0.2) &
#                         yield_trim < quantile(yield_trim, 0.6) & yield_trim > quantile(yield_trim, 0.2))
# dim(dat_tet_best)
# 
# dat_trim_best<- filter(dat_trim, yield_trim > quantile(yield_trim, 0.7) & 
#                         yield_kan < quantile(yield_kan, 0.6) & yield_kan > quantile(yield_kan, 0.2) &
#                         yield_tet < quantile(yield_tet, 0.6) & yield_tet > quantile(yield_tet, 0.2) &
#                         yield_cam < quantile(yield_cam, 0.6) & yield_cam > quantile(yield_cam, 0.2))
# dim(dat_trim_best)
# 
# # Filter strains for each antibiotic that is in top 50%, and then bottom 50% for the other three antibiotics
# dat_kan_best<- filter(dat_kan, yield_kan > quantile(yield_kan, 0.5) & 
#                         yield_cam < quantile(yield_cam, 0.5) &
#                         yield_tet < quantile(yield_tet, 0.5) &
#                         yield_trim < quantile(yield_trim, 0.5))
# 
# dat_tet_best<- filter(dat_tet, yield_tet > quantile(yield_tet, 0.5) & 
#                         yield_cam < quantile(yield_cam, 0.5) &
#                         yield_kan < quantile(yield_kan, 0.5) &
#                         yield_trim < quantile(yield_trim, 0.5))
# 
# dat_trim_best<- filter(dat_trim, yield_trim > quantile(yield_trim, 0.5) & 
#                         yield_cam < quantile(yield_cam, 0.5) &
#                         yield_kan < quantile(yield_kan, 0.5) &
#                         yield_tet < quantile(yield_tet, 0.5))
# 
# write.csv(dat_cam, file = "sorted_cam.csv")
# write.csv(dat_kan, file = "sorted_kan.csv")
# write.csv(dat_tet, file = "sorted_tet.csv")
# write.csv(dat_trim, file = "sorted_trim.csv")
# 
# # Combine the four dataframes into one df
# dat_all_best<- bind_rows(dat_cam_best, dat_kan_best, dat_tet_best, dat_trim_best)
# dim(dat_all_best)
# View(dat_all_best)
# write.csv(dat_all_best, file = "sorted_all.csv")
# 
# 
# # Also look at (i) normalised growth and (ii) genertaion time for the same strains. 
# # Create list with all strains in one list. 
# # Communicate this to Martin
# 
# # Check if strains with high yield are at the edges of plate!
# 
# #dat4<- dat3 %>% 
# #  arrange(desc(yield_cam)) %>% 
# #  arrange(desc(yield_kan)) %>% 
# #  arrange(desc(yield_tet)) %>% 
# #  arrange(desc(yield_trim)) 
# # View(dat4)
# 
# 
# 
