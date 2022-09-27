
1) Download all the files from the folder 'Mortality_analysis' into a folder in you computer and set it as the working directory for operation in R/RStudio.

    * List of countries: The file 'Countries_FSI.csv' contains the list of 39 countries considered in the mortality analysis.
    * Generation of Covariates: The covariate data used in the model was created using the script 'Covariate_generation.csv'. The output file 'X_centerscale.csv' is a 39x5 dataset whose columns are the standardized covariates for year 2013 in the follwoing order:
         
         i.  GDP year-on-year % change: download link: https://data.worldbank.org/indicator/NY.GDP.MKTP.KD.ZG .
         
         ii. Current healthcare expenditure % GDP: download link: https://data.worldbank.org/indicator/SH.XPD.CHEX.GD.ZS .
         
         iii.Carbon-dioxide emission in metric tonnes per capita: download link: https://data.worldbank.org/indicator/EN.ATM.CO2E.PC?name_desc=false.
         
         iv. Infant mortality per 1000 live births: download link: https://childmortality.org/data.
         
         v.  Human Development Index: download link: https://hdr.undp.org/data-center/documentation-and-downloads.
         
      The script 'Covariate_generation.R' uses as inputs the files 'Countries_FSI.csv' along with .csv files downloaded from the links above. The covariate data were downloaded on September 12, 2022.  
    
   * Generation of Response data: The age-at-death data used to create the pdfs and quantiles for 39 countries were downloaded from https://www.mortality.org/ on August 18, 2020. The codes in the script 'Response_generation_new.csv' file was used to construct the output files 'quant_all.csv' and 'density_all.csv' which are both 39x101 data frames whose rows represent the contries in the file 'Countries_FSI.csv'. The columns of 'quant_all.csv' represent the equidistant quantiles on the grid [0,1]. The columns of 'density_all.csv' represent the equidistant points on the support [20,110] for mortality distributions.


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
    
5) Run the codes in the script 'CV_folds_analysis.R'. To understand performance of the models better, we split the data into training/testing segments in 30 folds. The training split consists of 29 observation while testing split consists of 10 observations in each fold. The splits of 30 folds were picked randomly without replacement and stored in 'Folds_new.csv' file to be used for all models repeatedly. Then the respective Local Fréchet and Global Fréchet models were built on the training split and used for prediction on the testing set. The Mean Square Prediction Error is calculated on the testing set. The best bandwidths obtained by running the script 'Models_BW_Predictions.R' are used here. The inputs are:

    
    a) X_centerscale.csv
    
    b) quant_all.csv
    
    c) Folds_new.csv 
    
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

7) To run the computations for the table 5 in the paper run the codes in the script 'Table5_computation.R'. It sources the function 'Frechet_Rsquared.R'. The inputs are:

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
    

9)  To prepare a similar file as Folds_new.csv, run the codes in the script Folds_generate.csv. The test/train partition of folds will not be identical to
    the one used in our analysis. As output it produces:
    
    Folds_new.csv
