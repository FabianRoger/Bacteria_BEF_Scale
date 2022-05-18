


This repository contains the code for the simulations and analysis for the article by Gamfeldt et al. 

preliminary citation: 

> Scale and heterogeneity increase transgressive overyielding in biodiversity-ecosystem function experiments
Lars Gamfeldt, Fabian Roger, James Hagan, Martin Palm, Jonas Warringer, Anne Farewell

The analysis has three components: Data simulation, Experimental validation and Meta-analysis.

To reproduce the analysis for the data simulation and Fig. S6, run the script called:

+ 1_Data_simulation.Rmd

To reproduce the analysis for the experimental validation section, run the script called:

+ 2c_Experiment_analysis.nb.Rmd

You must run the script `2c_Experiment_analysis.nb.Rmd` twice, once with `Sens_analysis = FALSE` and once with `Sens_analysis = TRUE`, to produce the figures in the main text and the supplementary figures. Running this script will also produce Fig. S3 and S4.

In addition, there are two other scripts:

+ 2a_Join_and_export_raw_data.Rmd
+ 2b_Exploring_data_BSL2_6_7_ATCC_ABs_absolute.R

These scripts cannot be run as the necessary data are not provided. The first script documents how the raw data files from the OD reader in the experiment have been joined. The joined data are uploaded to figshare and downloaded in the Analysis.Rmd script. 

The second script documents how the strains have been selected to maximize variance in growth performance across environments. The script is for internal documentation. These data are not public. 

One the scripts above have been run, you can now reproduce Fig. 2, Fig. 3 and Fig. S7 using the script:

+ 3_Plot_Fig_2.Rmd

To reproduce the analysis for the meta-analysis, start with the following script:

+ 4a_Meta_analysis_data_assessment.R

This script reproduces the statistics regarding how many publications were searched, included etc. that are reported in the text.

To reproduce the meta-analysis results, run the following script:

+ 4b_Meta_analysis.R

This script reproduces all figures and tables associated with the meta-analysis from the raw data files, namely:

+ Fig. 4
+ Fig. S5

The script also reproduces additional statistics and values reported only in text.

In addition, 

All data for the meta-analysis are downloaded directly within the scripts from figshare.com. The two datasets and their associated meta-data can be found at: https://ndownloader.figshare.com/files/22918043



### Instruction to download a Github repo

#### with git

in the Terminal:

```cd path/to/local/folder``` 

(on your computer - the folder were you want the repository to live) command on Windows might differ. 


```git clone https://github.com/FabianRoger/Bacteria_BEF_Scale.git```

This should download the directory. 

#### without git
If you don't have git installed, you can download the repository as zip file and save it locally. 

--> Clone or download (green button top right)
--> Download Zip

then save and extract the zip where you want the director to live. 