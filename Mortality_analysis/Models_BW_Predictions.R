
## In this script we select bandwidth for each covariate and generate predictions
## using Local Frechet regression. Also we generate predictions from Global Frechet
## regression considering all covariates. 

## lastly we save the bandwidths and model predictions as densities/quantiles for
## further analysis.

# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

# read the covariate data:
X_ctr<- read.csv("X_centerscale.csv", header = T)
country_<- X_ctr[,1]
X_ctr <- X_ctr[,-1]
rownames(X_ctr)<- country_

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- quant_all[,1]
quant_all <- quant_all[,-1]

# define the number of quantiles
m<- ncol(quant_all)

# define the number of observations
n<- nrow(X_ctr)

# support for densities
dSup = seq(20, 110, length.out = m)

# support for quantiles
qSup = seq(0, 1, length.out = m)

###################################
## Local Frechet on each covariate 
###################################

library('frechet')
set.seed(1233)

# for gdp YoY % change
lf_gdpc<- LocDenReg(xin=as.vector(X_ctr[,"GDPC"]), qin=as.matrix(quant_all), 
                    optns = list(qSup = qSup, dSup=dSup,lower=20, upper=110))
h_gdpc<- lf_gdpc$optns$bwReg # storing the best bandwidth 
write.csv(h_gdpc, "LF_GDPC_BW.csv")

# save the predicted quantiles of mortality distributions
write.csv(lf_gdpc$qout, "LF_GDPC_Qpred.csv")
# save the predicted densities of mortality distributions
write.csv(lf_gdpc$dout, "LF_GDPC_Dpred.csv")

# for health care expenditure
lf_hce<- LocDenReg(xin=as.vector(X_ctr[,"HCE"]), qin=as.matrix(quant_all), 
                     optns= list(qSup = qSup, dSup=dSup, lower=20, upper=110))
h_hce<- lf_hce$optns$bwReg  # storing the best bandwidth 
write.csv(h_hce, "LF_HCE_BW.csv")

# save the predicted quantiles of mortality distributions
write.csv(lf_hce$qout, "LF_HCE_Qpred.csv")
# save the predicted densities of mortality distributions
write.csv(lf_hce$dout, "LF_HCE_Dpred.csv")

# for CO2 emissions in tonnes per capita
lf_co2e<- LocDenReg(xin=as.vector(X_ctr[,"CO2E"]), qin=as.matrix(quant_all), 
                    optns= list(qSup = qSup, dSup=dSup, lower=20, upper=110))
h_co2e<- lf_co2e$optns$bwReg  # storing the best bandwidth 
write.csv(h_co2e, "LF_CO2E_BW.csv")

# save the predicted quantiles of mortality distributions
write.csv(lf_co2e$qout, "LF_CO2E_Qpred.csv")
# save the predicted densities of mortality distributions
write.csv(lf_co2e$dout, "LF_CO2E_Dpred.csv")


# for Infant mortality per 1000 births 
lf_im<-LocDenReg(xin=as.vector(X_ctr[,"IM"]), qin=as.matrix(quant_all), 
                    optns=list(qSup=qSup, dSup=dSup,lower=20, upper=110))
h_im<- lf_im$optns$bwReg # storing the best bandwidth 
write.csv(h_im, "LF_IM_BW.csv")

# save the predicted quantiles of mortality distributions
write.csv(lf_im$qout, "LF_IM_Qpred.csv")

# save the predicted densities of mortality distributions
write.csv(lf_im$dout, "LF_IM_Dpred.csv")


# for HDI 
lf_hdi<-LocDenReg(xin=as.vector(X_ctr[,"HDI"]), qin=as.matrix(quant_all), 
                  optns=list(qSup=qSup, dSup=dSup,lower=20, upper=110))
h_hdi<- lf_hdi$optns$bwReg # storing the best bandwidth 
write.csv(h_hdi, "LF_HDI_BW.csv")

# save the predicted quantiles of mortality distributions
write.csv(lf_hdi$qout, "LF_HDI_Qpred.csv")

# save the predicted densities of mortality distributions
write.csv(lf_hdi$dout, "LF_HDI_Dpred.csv")


# creating Global Fréchet prediction
gf_pred<- GloDenReg(xin = as.matrix(X_ctr), qin= as.matrix(quant_all),
                  optns= list(qSup = qSup, dSup=dSup, lower=20, upper=110))

# save the predicted quantiles of mortality distributions
write.csv(gf_pred$qout, "GF_Qpred.csv")

# save the predicted densities of mortality distributions
write.csv(gf_pred$dout, "GF_Dpred.csv")

