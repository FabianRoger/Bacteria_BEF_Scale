
This repository contains the code for the simulations and analysis for the article by Gamfeldt et al. 

preliminary citation: 

> Scaling-up the biodiversity-ecosystem functioning relationship: The effect of environmental heterogeneity on transgressive overyielding (2022, Oikos)
Lars Gamfeldt, James G. Hagan, Martin Palm, Jonas Warringer, Anne Farewell, Fabian Roger

The analysis has three components: Data simulation, Experimental validation and Meta-analysis. All data are downloaded directly from within the scripts via figshare.

To reproduce the analysis for the data simulation and Fig. S6, run the script called:

+ 01_Data_simulation.Rmd

To reproduce the analysis for the experimental validation section, run the script called:

+ 02c_Experiment_analysis.nb.Rmd

This script will directly download the necessary data from a figshare repository which can be found at: https://doi.org/10.6084/m9.figshare.12279884.v1.

You must run the script `2c_Experiment_analysis.nb.Rmd` twice, once with `Sens_analysis = FALSE` and once with `Sens_analysis = TRUE`, to produce the figures in the main text and the supplementary figures. Running this script will also produce Fig. S3 and S4.

In addition, there are two other scripts:

+ 02a_Join_and_export_raw_data.Rmd
+ 02b_Exploring_data_BSL2_6_7_ATCC_ABs_absolute.R

These scripts cannot be run as the necessary data are not provided. The first script documents how the raw data files from the OD reader in the experiment have been joined. The joined data are uploaded to figshare and downloaded in the `02c_Experiment_analysis.nb.Rmd` script. 

The second script documents how the strains have been selected to maximize variance in growth performance across environments. The script is for internal documentation. These data are not public. 

Once the scripts above have been run, you can now reproduce Fig. 2, Fig. 3 and Fig. S7 using the script:

+ 03_Plot_Figs_2_3_S7.Rmd

To reproduce the analysis for the meta-analysis, start with the following script:

+ 04a_Meta_analysis_data_assessment.R

This script reproduces the statistics regarding how many publications were searched, included etc. that are reported in the text.

To reproduce the meta-analysis results, run the following script:

+ 04b_Meta_analysis.R

This script reproduces all figures and tables associated with the meta-analysis from the raw data files, namely:

+ Fig. 4
+ Fig. S5

The script also reproduces additional statistics and values reported only in text.

All data for the meta-analysis are downloaded directly within the scripts from figshare.com. The two datasets and their associated meta-data can be found at: https://ndownloader.figshare.com/files/22918043.
