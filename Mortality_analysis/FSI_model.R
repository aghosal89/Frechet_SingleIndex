
#################################
## Codes to Build FSI model
#################################

# set working directory~
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/FSI/Mortality_all")

# read the list of countries for our model~
country_<- read.csv("Countries_FSI.csv", header = T)[,2]

# read the covariate data:
X_ctr<- read.csv("X_centerscale.csv", header = T)
country_<- X_ctr[,1]
X_ctr <- X_ctr[,-1]
rownames(X_ctr)<- country_

# Read the response data as quantiles
quant_all <- read.csv("quant_all.csv", header = T)
rownames(quant_all)<- quant_all[,1]
quant_all <- quant_all[,-1]

# this function computes the euclidean norm of the vector x:
e_norm <- function(x) sqrt(sum(x^2))

# to find the bandwidth range for analysis~
h_max = max(apply(X_ctr, 1, e_norm))

## to find the lowest possible value for bandwidth h
metric_v_temp <- matrix(NA, 40 ,1) 
mv <-matrix(NA, 40 ,1)

for (j in 1:40) {
  for (i in 1:40) {
    if(i!=j)
      # computing the euclidean distance between rows the standardized X matrix
      mv[i]<- e_norm(X_ctr[j,] - X_ctr[i,]) 
  }
  metric_v_temp[j] = min(mv[-j]) # taking minimum distance 
}

h_min <- min(metric_v_temp)

# the sequence of bandwidths to optimize over ~
h = exp(seq(log(h_min),log(h_max), length.out = 10))

# define the length of quantiles
m<-ncol(quant_all)

# number of observations 
n<- nrow(X_ctr)

qSup<- seq(0,1, length.out = m)

library('numbers')
library('frechet')
library('fdadensity')

#######################################################
## Choosing bandwidth by leave-one-out Cross-Validation
#######################################################

mspe_l1ocv <- matrix(NA, length(h), 1) 
pe_temp <- matrix(0, n, 1)

for (k in 1:length(h)) {
  print(k)
  for (i in 1:n) {
    print(i)
    # select the quentiles in train and test splits of fold
    q_in <- quant_all[-i,] # training split
    q_out<- quant_all[i,]  # testing split
    
    # select the covariate observations in train and test splits of fold
    x_in <- X_ctr[-i,]  # training split
    x_out<- X_ctr[i,]   # testing split
    
    # fitting frechet single index model
    tempMatrix <- FSIDenReg(as.matrix(x_in), qSup, as.matrix(q_in), h[k], 
                            kern="gauss", Xout= as.matrix(x_out), 2)
    
    # computing MSPE on out-sample
    pe_temp[i]<- fdadensity:::trapzRcpp(X = qSup, Y = (as.vector(tempMatrix$Yout) - as.numeric(q_out))^2)
  }
  
  mspe_l1ocv[k,] <- mean(pe_temp)
}

h_fsi <- h[which.min(mspe_l1ocv)]

write.csv(h_fsi, 'FSI_bw.csv') 







###########################################
## Estimation of the single index parameter 
###########################################

nstart=5
tempMatrix <- FSIDenReg(as.matrix(X_ctr), qSup, as.matrix(quant_all), 
                        h_fsi, kern="gauss", Xout= as.matrix(X_ctr), nsp= nstart)

theta_hat <- tempMatrix$thetaHat

# theta_hat <- c(0.0774, 0.7484, 0.0327, 0.04754, 0.6562)
write.csv(theta_hat, "Theta_Hat.csv")

# generate the FSI model predictions
# obtain the frechet single index regression for the chosen bandwidth above
fsi_pred <- LocDenReg(xin = (as.matrix(X_ctr)%*% theta_hat), qin= as.matrix(quant_all),
                  optns= list(bwReg= h_fsi, qSup = qSup, dSup=dSup, lower=20, upper=110))

# save model predictions from Frechet Single Index model
write.csv(fsi_pred$qout, 'FSI_Qpred.csv')

# save the predicted densities of mortality distributions
write.csv(fsi_pred$dout, "FSI_Dpred.csv")


####################################
## Cross-Validation across 30 folds 
####################################

fold <- read.csv("Folds.csv", header = T)

mspe_kfcv <- matrix(NA, nrow(fold), 1)
pe_outfold   <- matrix(0, 30, 1)   # to store the Wn for testing set after theta optimization
pe_infold   <- matrix(0, 30, 1)   # to store minimized Wn for training set theta optimization
thetahat_fold <- matrix(0, 30, 5)  # store the theta estimate for each training split
nstart=5

for (j1 in 1:nrow(fold)) {
  print(j1)
  q_in <- quant_all[-as.numeric(fold[j1,]),]
  q_out<- quant_all[as.numeric(fold[j1,]),]
  
  x_in <- X_ctr[-as.numeric(fold[j1,]),]
  x_out<- X_ctr[as.numeric(fold[j1,]),]
  
  tempMatrix <- FSIDenReg(as.matrix(x_in), qSup, as.matrix(q_in), h=h_fsi, 
                          kern="gauss", Xout= as.matrix(x_out), nsp=nstart)
  thetahat_fold[j1, ] <- tempMatrix$thetaHat
  pe_infold[j1] <- tempMatrix$fnvalue
  
  pe_outfold[j1] <- mean(sapply(1:nrow(tempMatrix$Yout), 
                     function(i) fdadensity:::trapzRcpp(X= qSup, 
                      Y = (tempMatrix$Yout[i, ]- as.numeric(q_out[i, ]))^2)))
}

write.csv(pe_outfold, "FSI_MSPE_folds.csv")



