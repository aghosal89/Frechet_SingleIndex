
#########################################
# getting R-square oplus for the models:
#########################################

# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/LocalFrechetSpherical/Mortality_all")

library('frechet')

# reading the response data
Y <- read.csv("quant_all.csv", header = T)
Y<- as.matrix(Y[,-1])

# read covariate data
X<- read.csv("X_centerscale.csv", header = T)
X<- as.matrix(X[,-1])


qSup <- seq(0,1, length.out=101)

# For Global Frechet
frechet_Rsquared(Y=Y, X=X, tt=qSup, model="GF")*100

# For Local Frechet: HDI
h_hdi<- read.csv("LF_HDI_bw.csv", header = T)
h_hdi[,2]
frechet_Rsquared(Y=Y, X=(X[,"HDI"]), tt= qSup, h=h_hdi[,2], model="LF")*100


# For Local Frechet: Healthcare Expenditure as percentage of GDP
h_hcexp<- read.csv("LF_HCexp_bw.csv", header = T)
h_hcexp[,2]
frechet_Rsquared(Y=Y, X=(X[,"HCexp"]), tt= qSup, h=h_hcexp[,2], model="LF")*100


# For Local Frechet: GDP
h_gdp<- read.csv("LF_GDP_bw.csv", header = T)
h_gdp[,2]
frechet_Rsquared(Y=Y, X=(X[,"GDP"]), tt= qSup, h=h_gdp[,2], model="LF")*100


# For Local Frechet: Infant Mortality
h_infmt<- read.csv("LF_Infmt_bw.csv", header = T)
h_infmt[,2]
frechet_Rsquared(Y=Y, X=(X[,"Infmt"]), tt= qSup, h=h_infmt[,2], model="LF")*100


# For Local Frechet: CO2em
h_co2em<- read.csv("LF_CO2em_bw.csv", header = T)
h_co2em[,2]
frechet_Rsquared(Y=Y, X=(X[,"CO2em"]), tt= qSup, h=h_co2em[,2], model="LF")*100


# For Frechet Single Index
# read the bandwidth for model
h_fsi<- read.csv("FSI_bw.csv", header = T)
h_fsi[,2]

# read the theta hat for the model 
theta <- read.csv("Theta_Hat.csv", header = T)

frechet_Rsquared(Y=Y, X=(X%*%theta[,2]), tt= qSup, h=h_fsi[,2], model="LF")*100

