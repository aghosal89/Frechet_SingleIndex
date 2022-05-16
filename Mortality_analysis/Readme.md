
1) Download the all the files from the folder 'Mortality_analysis' into a folder in you computer and set it as the working directory for operation in R/RStudio.

2) Run the codes in the script 'Covariate_generation.R' to read and generate the covariate data for our model. 
This script takes as inputs the following files:
    
    a) Countries_FSI.csv:    contains list of all 40 countries in order appearing in the covariate data
    b) GDP_YoY_perc.csv:     contains time seeries data on GDP year-on-year % change for a superset of countries containing 40 countries in our model.
    c) HC_exp_percGDP.csv:   contains time seeries data on Healthcare expenditure % of GDP for a superset of countries containing 40 countries in our model.
    d) CO2_emissions_pc.csv: contains time seeries data on CO2 emissions per-capita for a superset of countries containing 40 countries in our model.
    e) infmort.csv:          contains time seeries data on Infant mortality per 1000 live births for a superset of countries containing 40 countries in our model.
    f) HumanDevelopmentIndex(HDI).csv: contains time seeries data on Human Development Index for a superset of countries containing 40 countries in our model.

As outputs it produces:

    a) X_centerscale.csv, the 40x5 design matrix for our model. 

3) Run the codes in the script 'Response_generation.R' to read and generate the response data in our model. This script takes as inputs the following files:
    
    3.1)    lt_australia.csv
    
    3.2)    lt_austria.csv
    
    3.3)    lt_belarus.csv  
    
    3.4)    lt_belgium.csv
    
    3.5)    lt_bulgaria.csv
    
    3.6)    lt_canada.csv
    
    3.7)    lt_chile.csv
    
    3.8)    lt_croatia.csv
    3.9)    lt_czech.csv
    3.10)   lt_denmark.csv
    3.11)   lt_estonia.csv 
    3.12)   lt_finland.csv
    3.13)   lt_france.csv
    3.14)   lt_germany.csv
    3.15)   lt_greece.csv
    3.16)   lt_hongkong.csv
    3.17)   lt_hungary.csv
    3.18)   lt_iceland.csv
    3.19)   lt_ireland.csv
    3.20)   lt_israel.csv
    3.21)   lt_italy.csv
    3.22)   lt_japan.csv
    3.23)   lt_latvia.csv
    3.24)   lt_lithuania.csv
    3.25)   lt_luxembourg.csv
    3.26)   lt_netherlands.csv
    3.27)   lt_newz.csv
    3.28)   lt_norway.csv
    3.29)   lt_poland.csv
    3.30)   lt_portugal.csv
    3.31)   lt_korea.csv
    3.32)   lt_russia.csv
    3.33)   lt_slovakia.csv
    3.34)   lt_slovenia.csv
    3.35)   lt_spain.csv
    3.36)   lt_sweden.csv
    3.37)   lt_switzerland.csv
    3.38)   lt_UK.csv
    3.39)   lt_USA.csv
    3.40)   lt_Ukraine.csv
    
    3.41)   Countries_FSI.csv

the files 3.1 - 3.40 above are life tables for the respective countries above for various years. Each contains columns for age and number of deaths 
in that age for that country. The file 'Countries_FSI.csv' is described above in 2). As outputs it generates the files:
    
    a) quant_all.csv  : a 40x101 matrix whose rows are the quantiles on an equidistant grid of length 101 on [0,1] of the mortality distribution of each country.
    b) density_all.csv : a 40x101 matrix whose rows are the mortality densities on an equidistant grid of age-range [20,110] for each country.


4) Run the codes in the script 'Models_BW_Predictions.R' to find the predicted quantiles/densities for all countries using Local Frechet regression using all 
covariates 'GDPC', 'HCE', 'CO2E', 'IM', 'HDI' and using Global Frechet regression. In the process we also find the best bandwidth choices for each Local Frechet 
regression models. The inputs for this script are:

    a) X_centerscale.csv
    b) quant_all.csv

The outputs are:

    a) LF_GDPC_BW.csv
    b) LF_GDPC_Qpred.csv
    c) LF_GDPC_Dpred.csv
    
    d) LF_HCE_BW.csv
    e) LF_HCE_Qpred.csv
    f) LF_HCE_Dpred.csv
    
    g) LF_CO2E_BW.csv
    h) LF_CO2E_Qpred.csv
    i) LF_CO2E_Dpred.csv
    
    j) LF_IM_BW.csv
    k) LF_IM_Qpred.csv
    l) LF_IM_Dpred.csv
    
    m) LF_HDI_BW.csv
    n) LF_HDI_Qpred.csv
    o) LF_HDI_Dpred.csv
    
    p) GF_Qpred.csv
    q) GF_Dpred.csv    
    
5) Run the codes in the script 'CV_folds_analysis.R'. Here we consider 30 folds of training/testing splits in data. In each fold we picked randomly 
without replacement 1/4 observations to be in testing split while rest 3/4 observations to be in the training split. Then the respective Local Frechet
or Global Frechet models are built on the training split and used for predicting the response on the testing set. The Mean Square Prediction Error is
measured on the testing set. The best bandwidths obtained by running the script Models_BW_Predictions.R are used here. The inputs are:

    a) X_centerscale.csv
    b) quant_all.csv
    c) Folds.csv 
    
    d) LF_GDPC_BW.csv
    e) LF_HCE_BW.csv
    f) LF_CO2E_BW.csv
    g) LF_IM_BW.csv
    h) LF_HDI_BW.csv
    
To run the codes source the function from the script 'LocWassRegAMP.R'.

The outputs are:

    a) GF_folds.csv
    b) LF_GDPC_folds.csv
    c) LF_HCE_folds.csv
    d) LF_CO2E_folds.csv
    e) LF_IM_folds.csv
    f) LF_HDI_folds.csv

6) To run the FSI model, run the codes in the script 'FSI_model.R', first we read the functions from the following scripts in your working directory:

    a) FSIAuxFunctions.R
    b) FSIDenReg.R
    c) LocWassRegAMP.R

The FSI_model.R script takes as input the following files: 

7) To run the computations for the table 5 in the paper run the codes in the script 'Table5_computation.R'. Need to source the function 'frechet_Rsquared.R'.

8) To generate the figures 7 - 10 in the paper run the codes in the script 'Plots.R'. 
