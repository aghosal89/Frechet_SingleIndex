

#######################################
# Computations for Table 5
# get Frechet R-squared for the models
#######################################

# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

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

# function to produce the R_oplus_squared for Wasserstein space

# INPUTS 1) Y    : an n by m matrix where each row describes the quantiles of the distribution. 
#        2) X    : an n by p vector of covariates for Global Frechet model, n by 1 vector for Local Frechet model.
#        3) tt   : an equidistant grid of length m on [0,1] on which the quantiles are computed.
#        4) h    : Bandwidth for the model, defaults to NA incase of Global Frechet model.
#        5) model: To specify "GF" or "LF" for Global Frechet or Local Frechet respectively.


frechet_Rsquared <- function(Y=NULL, X=NULL, h= NULL, tt=NULL, model=NULL) {
  
  if(model=="LF" & is.null(h)==TRUE) {
    error("For Local Frechet model a bandwidth required.")
  }
  
  if(model=="GF") {
    
    fr<- frechet:::GloWassReg(xin=(X), qin=(Y), xout = X, 
                              optns =list(qSup = tt, dSup= seq(20,110, length.out=101), lower=20, upper=110))
    
    qmean <- apply(fr$qout, 2, mean)  
    
    temp_n <- matrix(NA, dim(X)[1], 1)
    temp_d <- temp_n
    
    for (i in 1:dim(X)[1]) {
      temp_n[i]<- fdadensity:::trapzRcpp(X= tt, Y = (Y[i,]- fr$qout[i,])^2)
      temp_d[i]<- fdadensity:::trapzRcpp(X= tt, Y = (Y[i,]-qmean)^2)
    }
    
    rsq_oplus <- 1- (sum(temp_n)/sum(temp_d))
    return(rsq_oplus)
    
  } else if(model=="LF") {
    fr<- frechet:::LocWassReg(xin=as.vector(X), qin=Y, xout = as.vector(X), 
                              optns =list(bwReg=h, qSup = tt, dSup= seq(20,110, length.out=101), lower=20, upper=110))
    qmean <- apply(fr, 2, mean)  
    
    temp_n <- matrix(NA, length(X), 1)
    temp_d <- temp_n
    
    for (i in 1:length(X)) {
      temp_n[i]<- fdadensity:::trapzRcpp(X= tt, Y = (Y[i,]- fr[i,])^2)
      temp_d[i]<- fdadensity:::trapzRcpp(X= tt, Y = (Y[i,]- qmean)^2)
    }
    rsq_oplus <- 1- (sum(temp_n)/sum(temp_d))
    return(rsq_oplus)
  }
  
}

library("frechet")
library("fdadensity")

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
write.csv( df_mspe_fold,"MSPE_folds.csv")




