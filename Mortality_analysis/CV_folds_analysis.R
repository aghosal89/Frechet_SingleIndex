############
# Script #4
############

## This script divides the dataset into various folds then
## computes MSPE across folds for all models with Local Frechet regression with the 
## single covariates, Global Frechet regression with all covariates. 

# The performance of the FSI model is described in the script 'FSI_model.R'.

# Aritra's working directory~
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

# we read the folds for test-train splits for response data
fold <- read.csv("Folds.csv", header = T)

# define the number of quantiles
m<-ncol(quant_all)

# define the number of observations
n<- nrow(X_ctr)

# grid of support for densities 
dSup = seq(20, 110, length.out = m)

# grid of support for quantiles
qSup = seq(0, 1, length.out = m)

###################################
## Global Frechet on all covariates
###################################

library('frechet')

# creating matrix to store the mspe across folds
gf_folds <- matrix(NA, 30, 1) 
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),])
  fr <- frechet:::GloWassReg(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list( qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr$qout[s,]- q_out[s,])^2) 
  }
  
  gf_folds[i,] <- mean(pe)
}

write.csv(gf_folds, "GF_folds.csv")


# GDP Year-on-Year
# read the bandwidth
h_gdpc<- read.csv("LF_GDPC_BW.csv", header = T)[,2]

lf_gdpc_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),])
  fr <- LocWassRegAMP(xin = as.matrix(x_in[,"GDPC"]), qin= q_in, xout= x_out[,"GDPC"], 
                  optns= list(bwReg=h_gdpc, qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(as.numeric(fold[i,])))
  for(s in 1:length(fold[i,])) {
    pe[s]<- fdadensity:::trapzRcpp(X= qSup, Y=( fr[s,] - q_out[s,])^2)
  }
  lf_gdpc_folds[i,] <- mean(pe)
}
write.csv(lf_gdpc_folds, 'LF_GDPC_folds.csv')


# HDI
# Read the bandwidth
h_hdi<- read.csv("LF_HDI_BW.csv", header = TRUE)[,2]

lf_hdi_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),"HDI"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"HDI"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_hdi, qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,]-q_out[s,])^2)
  }
  lf_hdi_folds[i,] <- mean(pe)
}
write.csv(lf_hdi_folds, 'LF_HDI_folds.csv')


# Healthcare expenditure as percentage of GDP
h_hce <- read.csv("LF_HCE_BW.csv", header =  T)[,2]
lf_hce_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])), "HCE"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"HCE"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_hce, qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(as.numeric(fold[i,])))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  
  lf_hce_folds[i,] <- mean(pe)
}
write.csv(lf_hcexp_folds, 'LF_HCE_folds.csv')


# CO2 emissions in tonnes per capita
h_co2e <- read.csv("LF_CO2E_BW.csv", header = T)[,2]
lf_co2e_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),"CO2E"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"CO2E"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_co2e, qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  lf_co2e_folds[i,] <- mean(pe)
}
write.csv(lf_co2e_folds, "LF_CO2E_folds.csv")


# infant mortality per 1000 births
h_im<- read.csv("LF_IM_BW.csv", header=T)[,2]
lf_im_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])), "IM"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])), "IM"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_im, qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  
  lf_im_folds[i,] <- mean(pe)
}

write.csv(lf_im_folds, "LF_IM_folds.csv")



