
1) Download all the files from the folder 'Mortality_analysis' into a folder in you computer and set it as the working directory for operation in R/RStudio.


2) Run the codes in the script 'Covariate_generation.R' to read and process the covariate data for our model. 
This script takes as inputs the following files:
    
    
    a) Countries_FSI.csv:  contains list of all 39 countries in order appearing in the covariate data.
    
    
    b) GDP_YoY_perc.csv:  contains time seeries data on GDP year-on-year % change for a superset of countries containing 39 countries in our model.
    
    
    c) HC_exp_percGDP.csv:  contains time seeries data on Healthcare expenditure % of GDP for a superset of countries containing 39 countries in our model.
    
    
    d) CO2_emissions_pc.csv: contains time series data on CO2 emissions per-capita for a superset of countries containing 39 countries in our model.
    
    
    e) infmort.csv:   contains time series data on Infant mortality per 1000 live births for a superset of countries containing 39 countries in our model.
    
    
    f) HumanDevelopmentIndex(HDI).csv: contains time series data on Human Development Index for a superset of countries containing 39 countries in our model.
    

Sources:
The data on GDP year-on-year % change, Healthcare expenditure percentage of GDP, CO2 emissions in tonnes per capita, Infant mortality per 1000 live births are obtained from world bank website (https://data.worldbank.org/). The data on Human Development Index (HDI) is obtained from (https://hdr.undp.org/en/data). The data and the list of countries considered in our mortality analysis are obtained from the website 'mortality.org' (https://www.mortality.org/).

As output it produces:

    a) X_centerscale.csv, the 39x5 design matrix for our model whose columns are the covariates scaled and centered with each row/observation representing a country. We considered mortality analysis for the year 2013. In our analysis, we refer to these covariates as GDPC, HCE, CO2E, IM, HDI respectively in the order appearing above.


3) Run the codes in the script 'Response_generation.R' to read and process the response data in our model. This script takes as inputs the following files:
    
    
    1.   lt_australia.csv
    2.   lt_austria.csv
    3.   lt_belarus.csv  
    4.   lt_belgium.csv
    5.   lt_bulgaria.csv
    6.   lt_canada.csv
    7.   lt_chile.csv
    8.   lt_croatia.csv
    9.   lt_czech.csv
    10.  lt_denmark.csv
    11.  lt_estonia.csv 
    12.  lt_finland.csv
    13.  lt_france.csv
    14.  lt_germany.csv
    15.  lt_greece.csv
    16.  lt_hungary.csv
    17.  lt_iceland.csv
    18.  lt_ireland.csv
    19.  lt_israel.csv
    20.  lt_italy.csv
    21.  lt_japan.csv
    22.  lt_latvia.csv
    23.  lt_lithuania.csv
    24.  lt_luxembourg.csv
    25.  lt_netherlands.csv
    26.  lt_newz.csv
    27.  lt_norway.csv
    28.  lt_poland.csv
    29.  lt_portugal.csv
    30.  lt_korea.csv
    31.  lt_russia.csv
    32.  lt_slovakia.csv
    33.  lt_slovenia.csv
    34.  lt_spain.csv
    35.  lt_sweden.csv
    36.  lt_switzerland.csv
    37.  lt_UK.csv
    38.  lt_USA.csv
    39.  lt_Ukraine.csv  
    
    and finally,
    
    40.   Countries_FSI.csv

the files 1 - 39 above are life tables for the respective countries above for various years. The data obtained from (https://www.mortality.org/). Each contains columns for age and number of deaths in that age for that country. The file 'Countries_FSI.csv' is described above in 2). 

As outputs it generates the files:
    
    a) quant_all.csv  :  a 39x101 matrix whose rows are the quantiles on an equidistant grid of length 101 on [0,1] of the mortality distribution of each country.
    b) density_all.csv : a 39x101 matrix whose rows are the mortality densities on an equidistant grid of age-range [20,110] of length 101 for each country.

4) Run the codes in the script 'Models_BW_Predictions.R' to find the predicted quantiles/densities for all countries using Local Frechet regression using each of the covariates 'GDPC', 'HCE', 'CO2E', 'IM', 'HDI'; and also running Global Fréchet regression using all of these 5 covariates. In the process we also find the best bandwidth choices for each Local Frechet regression model (i.e. for Local Fréchet regression with each of the covariates). The inputs for this script are:

    a) X_centerscale.csv
    b) quant_all.csv

The outputs are:

    a) LF_GDPC_BW.csv    : contains the chosen bandwidth for the covariate GDPC.
    b) LF_GDPC_Qpred.csv : contains the Local Frechet regression predicted quantiles for 2013 for all countries using GDPC as the covariate.
    c) LF_GDPC_Dpred.csv : contains the Local Frechet regression predicted densities for 2013 for all countries using GDPC as the covariate.
    
    d) LF_HCE_BW.csv    : contains the chosen bandwidth for the covariate HCE.
    e) LF_HCE_Qpred.csv : contains the Local Frechet regression predicted quantiles for 2013 for all countries using HCE as the covariate.
    f) LF_HCE_Dpred.csv : contains the Local Frechet regression predicted densities for 2013 for all countries using HCE as the covariate.
    
    g) LF_CO2E_BW.csv    : contains the chosen bandwidth for the covariate CO2E.
    h) LF_CO2E_Qpred.csv : contains the Local Frechet regression predicted quantiles for 2013 for all countries using CO2E as the covariate.
    i) LF_CO2E_Dpred.csv : contains the Local Frechet regression predicted densities for 2013 for all countries using CO2E as the covariate.
    
    j) LF_IM_BW.csv    : contains the chosen bandwidth for the covariate IM.
    k) LF_IM_Qpred.csv : contains the Local Frechet regression predicted quantiles for 2013 for all countries using IM as the covariate.
    l) LF_IM_Dpred.csv : contains the Local Frechet regression predicted densities for 2013 for all countries using IM as the covariate.
    
    m) LF_HDI_BW.csv    : contains the chosen bandwidth for the covariate HDI.
    n) LF_HDI_Qpred.csv : contains the Local Frechet regression predicted quantiles for 2013 for all countries using HDI as the covariate.
    o) LF_HDI_Dpred.csv : contains the Local Frechet regression predicted densities for 2013 for all countries using HDI as the covariate.
    
    p) GF_Qpred.csv  : contains the Global Frechet regression predicted quantiles for 2013 for all countries using all the covariates.
    q) GF_Dpred.csv  : contains the Global Frechet regression predicted densities for 2013 for all countries using all the covariates.
    
5) Run the codes in the script 'CV_folds_analysis.R'. To understand performance of the models better, we split the data into training/testing segments in 30 folds. The training split consists of 30 observation while testing split consists of 10 observations in each fold. The splits of 30 folds were picked randomly without replacement and stored in 'Folds.csv' file to be used for all models repeatedly. Then the respective Local Frechet and Global Frechet models were built on the training split and used for prediction on the testing set. The Mean Square Prediction Error is calculated on the testing set. The best bandwidths obtained by running the script 'Models_BW_Predictions.R' are used here. The inputs are:

    a) X_centerscale.csv
    b) quant_all.csv
    c) Folds.csv 
    
    d) LF_GDPC_BW.csv
    e) LF_HCE_BW.csv
    f) LF_CO2E_BW.csv
    g) LF_IM_BW.csv
    h) LF_HDI_BW.csv
    
To run the codes source the function from the script 'LocWassRegAMP.R'.

The outputs: each file below contains 30x1 matrix whose elements are the Mean Square Prediction Error for the models in the 30 folds.

    a) GF_folds.csv       : MSPE for Global Frechet using all covariates in 30 folds.
    
    b) LF_GDPC_folds.csv  : MSPE for Local Frechet using GDPC in 30 folds.
    
    c) LF_HCE_folds.csv   : MSPE for Local Frechet using HCE in 30 folds.
    
    d) LF_CO2E_folds.csv  : MSPE for Local Frechet using CO2E in 30 folds.
    
    e) LF_IM_folds.csv    : MSPE for Local Frechet using IM in 30 folds.
    
    f) LF_HDI_folds.csv   : MSPE for Local Frechet using HDI in 30 folds.

6) To run the FSI model, find the best bandwidth, estimate of the index parameter, obtain the Mean Square Prediction Error in 30 folds, run the codes in the script 'FSI_model.R', which works by sourcing the functions from the following scripts in your working directory:

    a) FSIAuxFunctions.R
    
    b) FSIDenReg.R
    
    c) LocWassRegAMP.R

The FSI_model.R script takes as inputs the following files: 

    1. X_centerscale.csv
    
    2. quant_all.csv
    
    3. Countries_FSI.csv
    
The outputs:

    1. FSI_bw.csv       : the best chosen bandwidth for FSI model.
    
    2. Theta_Hat.csv    : Estimate of the index parameter using the best bandwidth above.
    
    3. FSI_Qpred.csv    : Predicted quantiles by FSI model by best bandwidth and estimated index parameter above. 
    
    4. FSI_Dpred.csv    : Predicted densities by FSI model by best bandwidth and estimated index parameter above.
    
    5. FSI_MSPE_folds.csv : Mean Square Prediction Error in 30 folds by the FSI model using the best bandwidth above. 

7) To run the computations for the table 5 in the paper run the codes in the script 'Table5_computation.R'. It sources the function 'frechet_Rsquared.R'. The inputs are:

    1. Countries_FSI.csv
    2. X_centerscale.csv
    3. quant_all.csv
    
    4. LF_HDI_BW.csv
    5. LF_HCE_BW.csv
    6. LF_GDPC_BW.csv
    7. LF_IM_BW.csv
    8. LF_CO2E_BW.csv
    9. FSI_bw.csv
    10. Theta_Hat.csv
    
    11. GF_folds.csv
    12. LF_HDI_folds.csv
    13. LF_HCE_folds.csv
    14. LF_GDPC_folds.csv
    15. LF_IM_folds.csv
    16. LF_CO2E_folds.csv

    17. FSI_MSPE_folds.csv

Outputs: This script does not produce outputs to be saved for further use. It is used to compute the numbers in table 5 of the paper. 

8) To generate the figures 7 - 10 in the paper run the codes in the script 'Plots.R'. The following are the inputs:

    a) density_all.csv
    
    b) GF_Dpred.csv
    
    c) LF_HDI_Dpred.csv
    
    d) FSI_Dpred.csv
    
    e) MSPE_folds.csv
    
    f) LF_HDI_folds.csv
    
    g) LF_HCE_folds.csv
    
    h) FSI_MSPE_folds.csv
    
    i) GF_folds.csv
    
    j) X_centerscale.csv
    
    j) quant_all.csv
    
    k) Theta_Hat.csv
    
    l) FSI_bw.csv
    

9)  To prepare a similar file as Folds.csv, run the codes in the script Folds_generate.csv. The test/train partition of folds will not be identical to
    the one used in our analysis. As output it produces:
    
    Folds_new.csv
