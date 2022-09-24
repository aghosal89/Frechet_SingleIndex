
#######################################
# Computations for Table 5
# get Frechet R-squared for the models
#######################################

# set working directory~

# files to be sourced~
source("frechet_Rsquared.R")

# read libraries~
library("frechet")
library("fdadensity")

# read the list of countries for our model~
country_<- read.csv("Countries_FSI.csv", header = T)[,2]

# read the covariate data:
X_ctr<- read.csv("X_centerscale.csv", header = T)
X_ctr <- as.matrix(X_ctr[,-1])
rownames(X_ctr)<- country_

# number of observations
n <- nrow(X_ctr)

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- country_
quant_all <- as.matrix(quant_all[,-1])

# length of quantiles 
m <- ncol(quant_all)

# support for quantiles
qSup <- seq(0,1, length.out=m)

# support for densities
dSup <- seq(20, 110, length.out=m)

# For Global Frechet
frechet_Rsquared(Y=quant_all, X=X_ctr, tt=qSup, model="GF")*100

# For Local Frechet: HDI
h_hdi<- read.csv("LF_HDI_BW.csv", header = T)[,2]
frechet_Rsquared(Y=quant_all, X=(X_ctr[,"HDI"]), tt= qSup, h=h_hdi, model="LF")*100

# For Local Frechet: Healthcare Expenditure as percentage of GDP
h_hce<- read.csv("LF_HCE_BW.csv", header = T)[,2]
frechet_Rsquared(Y=quant_all, X=(X_ctr[,"HCE"]), tt= qSup, h=h_hce, model="LF")*100

# For Local Frechet: GDP
h_gdpc<- read.csv("LF_GDPC_BW.csv", header = T)[,2]
frechet_Rsquared(Y=quant_all, X=(X_ctr[,"GDPC"]), tt= qSup, h=h_gdpc, model="LF")*100

# For Local Frechet: Infant Mortality
h_im<- read.csv("LF_IM_BW.csv", header = T)[,2]
frechet_Rsquared(Y=quant_all, X=(X_ctr[,"IM"]), tt= qSup, h=h_im, model="LF")*100


# For Local Frechet: CO2em
h_co2e<- read.csv("LF_CO2E_BW.csv", header = T)[,2]
frechet_Rsquared(Y=quant_all, X=(X_ctr[,"CO2E"]), tt= qSup, h=h_co2e, model="LF")*100


# For Frechet Single Index
# read the bandwidth for model
h_fsi<- read.csv("FSI_bw.csv", header = T)[,2]

# read the theta hat for the model 
theta <- read.csv("Theta_Hat.csv", header = T)[,2]

frechet_Rsquared(Y=quant_all, X=(X_ctr%*%theta), tt= qSup, h=h_fsi, model="LF")*100


######################################################################
## Computing the mean and standard deviation of the MSPE across folds 
######################################################################

# creating and saving dataset for MSPE variation across folds~

gf_folds<- read.csv("GF_folds.csv", header = T)
lf_hdi_folds <- read.csv("LF_HDI_folds.csv", header= T)
lf_hce_folds<- read.csv("LF_HCE_folds.csv", header =T)
lf_gdpc_folds <- read.csv("LF_GDPC_folds.csv",header=T)
lf_im_folds <- read.csv("LF_IM_folds.csv", header = T)
lf_co2e_folds <- read.csv("LF_CO2E_folds.csv", header = T)
fsi_folds <- read.csv("FSI_MSPE_folds.csv", header = T)

# creating and saving dataset for MSPE variation across folds~

df_mspe_folds <- rbind(data.frame(Model= 'GF',MSPE =gf_folds[,2]), data.frame(Model='LF(GDPC)', MSPE=lf_gdpc_folds[,2]),
                      data.frame(Model= 'LF(HCE)',MSPE=lf_hce_folds[,2]),
                      data.frame(Model= 'LF(CO2E)', MSPE=lf_co2e_folds[,2]), 
                      data.frame(Model= 'LF(IM)', MSPE=lf_im_folds[,2]), 
                      data.frame(Model= 'LF(HDI)', MSPE=lf_hdi_folds[,2]), 
                      data.frame(Model= 'FSI', MSPE = fsi_folds[,2]))

# getting mean and standard deviation of MSPE across folds 
# Global Frechet
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "GF")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "GF")[,2])


# Local Frechet: Human Development Index
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(HDI)")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(HDI)")[,2])


# Local Frechet: Health care expenditure as percentage of GDP
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(HCE)")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(HCE)")[,2])


# Local Frechet: GDP Year on year % change
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(GDPC)")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(GDPC)")[,2])


# Local Frechet: Infant Mortality in 1000 birth live births
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(IM)")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(IM)")[,2])


# Local Frechet: CO2 emissions in tonnes per capita
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(CO2E)")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "LF(CO2E)")[,2])


# Frechet Single Index model
mean(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "FSI")[,2])
sd(subset(df_mspe_folds, df_mspe_folds[,'Model'] == "FSI")[,2])


# save the MSPE distributions across folds for the models
write.csv( df_mspe_folds,"MSPE_folds.csv")




