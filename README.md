


This repository contains the code for the simulations and analysis for the article by Gamfeldt et al. 

preliminary citation: 

> Biodiversity and ecosystem functioning across gradients in spatial scale
Lars Gamfeldt, Fabian Roger, Martin Palm, James Hagan, Jonas Warringer, Anne Farewell


To reproduce the analysis, download the repository and execute the scripts 

+ Simulations.Rmd
+ Analysis.Rmd
+ Multipanel_Fig_2.Rmd

Th repository contains two further scripts: 

+ Join_and_export_raw_data.Rmd
+ exploring data BSL2_6_7_ATCC_ABs_absolute.R


These scripts cannot be run as the necessary data are not provided. The first script documents how the raw data files from the OD reader in the experiment have been joined. The joined data are uploaded to figshare and downloaded in the Analysis.Rmd script. 

The second script documents how the strains have been selected to maximize variance in growth performance across environments. The script is for internal documentation. These data are not public. 

### Instruction to download a Github repo

#### with git

in the Terminal:

```cd path/to/local/folder``` 

(on your computer - the folder were you want the repository to live) command on Windosw might differ. 


```git clone https://github.com/FabianRoger/Bacteria_BEF_Scale.git```

This should download the directory. 

#### without git
If you don't have git installed, you can download the repository as zip file and save it locally. 

--> Clone or download (green button top right)
--> Download Zip

then save and extract the zip where you want the director to live. 