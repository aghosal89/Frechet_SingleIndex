
# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/LocalFrechetSpherical/Mortality_all")

# read the covariate data:
X_ctr<- read.csv("X_centerscale.csv", header = T)
country_<- X_ctr[,1]
X_ctr <- X_ctr[,-1]
rownames(X_ctr)<- country_

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- quant_all[,1]
quant_all <- quant_all[,-1]

##### Computation of MSPE across folds for all models except FSI~
# we read the folds for test-train splits for response data
fold <- read.csv("Folds.csv", header = T)

dSup = seq(20, 110, length.out = 101)
qSup = seq(0, 1, length.out = 101)

######################################
## Global Frechet on all covariates
######################################

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

write.csv(gf_folds, "GF_mspe_folds.csv")

###################################
## Local Frechet on each covariate 
###################################

set.seed(1233)
# for gdp YoY % change
lf_gdp<- LocDenReg(xin=as.vector(X_ctr[,"GDP"]), xout=as.vector(X_ctr[,"GDP"]),
              qin=as.matrix(quant_all), optns = list(qSup = qSup,lower=20, upper=110))
h_gdp<- lf_gdp$optns$bwReg # storing the best bandwidth 
write.csv(h_gdp, "LF_GDP_BW.csv")

# for health care expenditure
lf_hcexp<- LocDenReg(xin=as.vector(X_ctr[,"HCexp"]), qin=as.matrix(quant_all), 
                     optns= list(qSup = qSup,lower=20, upper=110))
h_hcexp<- lf_hcexp$optns$bwReg # storing the best bandwidth 
write.csv(h_hcexp, "LF_HCEXP_BW.csv")

# for CO2 emissions in tonnes per capita
lf_co2em<-LocDenReg(xin=as.vector(X_ctr[,"CO2em"]), qin=as.matrix(quant_all), 
                    optns= list(qSup = qSup, lower=20, upper=110))
h_co2em<- lf_co2em$optns$bwReg  # storing the best bandwidth 
write.csv(h_co2em, "LF_CO2em_BW.csv")

# for Infant mortality per 1000 births 
lf_infmt<-LocDenReg(xin=as.vector(X_ctr[,"Infmt"]), qin=as.matrix(quant_all), 
                      optns=list(qSup=qSup,lower=20, upper=110))
h_infmt<- lf_infmt$optns$bwReg # storing the best bandwidth 
write.csv(h_infmt, "LF_Infmt_BW.csv")

# for HDI 
lf_hdi<-LocDenReg(xin=as.vector(X_ctr[,"HDI"]), qin=as.matrix(quant_all), 
                  optns=list(qSup=qSup,lower=20, upper=110))
h_hdi<- lf_hdi$optns$bwReg # storing the best bandwidth 
write.csv(h_hdi, "LF_HDI_BW.csv")

# GDP Year-on-Year
# read the bandwidth
h_gdp<- read.csv("LF_GDP_bw.csv", header = T)

lf_gdp_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),])
  fr <- LocWassRegAMP(xin = as.matrix(x_in[,"GDP"]), qin= q_in, xout= x_out[,"GDP"], 
                  optns= list(bwReg=h_gdp[,2], qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(as.numeric(fold[i,])))
  for(s in 1:length(fold[i,])) {
    pe[s]<- fdadensity:::trapzRcpp(X= qSup, Y=( fr[s,] - q_out[s,])^2)
  }
  lf_gdp_folds[i,] <- mean(pe)
}
write.csv(lf_gdp_folds, 'LF_GDP_folds.csv')


# HDI
# Read the bandwidth
h_hdi<- read.csv("LF_HDI_bw.csv", header = TRUE)

lf_hdi_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),"HDI"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"HDI"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_hdi[,2], qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,]-q_out[s,])^2)
  }
  lf_hdi_folds[i,] <- mean(pe)
}
write.csv(lf_hdi_folds, 'LF_HDI_folds.csv')


# Healthcare expenditure as percentage of GDP
h_hcexp <- read.csv("LF_HCExp_bw.csv", header =  T)
lf_hcexp_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])), "HCexp"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"HCexp"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_hcexp[,2], qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(as.numeric(fold[i,])))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  
  lf_hcexp_folds[i,] <- mean(pe)
}
write.csv(lf_hcexp_folds, 'LF_hcexp_folds.csv')


# CO2 emissions in tonnes per capita
h_co2em <- read.csv("LF_CO2em_bw.csv", header = T)
lf_co2em_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])),"CO2em"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])),"CO2em"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_co2em[,2], qSup = qSup, dSup=dSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  
  lf_co2em_folds[i,] <- mean(pe)
}
write.csv(lf_co2em_folds, "LF_co2em_folds.csv")


# infant moertality per 1000 births
h_infmt<- read.csv("LF_Infmt_bw.csv", header=T)
lf_infmt_folds <- matrix(NA, 30, 1)
for (i in 1:nrow(fold)) {
  
  q_in <- as.matrix(quant_all[-as.numeric(fold[i,]), ])
  q_out <- as.matrix(quant_all[as.numeric(fold[i,]), ])
  
  x_in <- as.matrix(X_ctr[-as.numeric(as.numeric(fold[i,])), "Infmt"])
  x_out <- as.matrix(X_ctr[as.numeric(as.numeric(fold[i,])), "Infmt"])
  fr <- LocWassRegAMP(xin = as.matrix(x_in), qin= q_in, xout= x_out, 
                  optns= list(bwReg=h_infmt[,2], qSup = qSup, lower=20, upper=110))
  
  pe <- matrix(NA, length(fold[i,]))
  for(s in 1:length(fold[i,])) {
    pe[s] <- fdadensity:::trapzRcpp(X= qSup, Y= (fr[s,] - q_out[s,])^2)
  }
  
  lf_infmt_folds[i,] <- mean(pe)
}

write.csv(lf_infmt_folds, "LF_infmt_folds.csv")









