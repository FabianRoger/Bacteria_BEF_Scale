


This repository contains the code for the simulations and analysis for the article by Gamfeldt et al. 

preliminary citation: 

> Biodiversity and ecosystem functioning across gradients in spatial scale
Lars Gamfeldt, Fabian Roger, Martin Palm, James Hagan, Jonas Warringer, Anne Farewell


To reproduce the analysis for the simulations and the bacteria microcosm experiment, download the repository and execute the scripts:

+ Simulations.Rmd
+ Analysis.Rmd
+ Multipanel_Fig_2.Rmd

Th repository contains two further scripts: 

+ Join_and_export_raw_data.Rmd
+ exploring data BSL2_6_7_ATCC_ABs_absolute.R


These scripts cannot be run as the necessary data are not provided. The first script documents how the raw data files from the OD reader in the experiment have been joined. The joined data are uploaded to figshare and downloaded in the Analysis.Rmd script. 

The second script documents how the strains have been selected to maximize variance in growth performance across environments. The script is for internal documentation. These data are not public. 


To reproduce the analysis for the meta-analysis, download the repository and execute the script:

+ Bef_scale_synthesis_script.R

This script reproduces all figures and tables associated with the meta-analysis from the raw data files, namely:

+ Fig. 3Ai - iii
+ Fig. 3B
+ Fig. S8
+ Fig. S9
+ Fig. S10
+ Table S1
+ Table S2
+ Table S3
+ Table S4

The script also reproduces additional statistics and values reported only in text

In addition, the script examines the data assessment file which was used to select papers for the meta-analysis via various inclusion/exclusion criteria.

All data are downloaded directly within the script from figshare.com. The meta-data for the two datasets used can be found at: https://ndownloader.figshare.com/files/22918043



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