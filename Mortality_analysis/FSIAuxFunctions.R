## This file contains auxiliary functions for use with FSIDenReg.R

polar2cart <- function(eta, r  = 1) {
  s <- length(eta)
  theta <- vector(length = s+1)
  
  if(s==4) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*cos(eta[4])
    theta[2]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])*sin(eta[4])
    theta[3]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
    theta[4]<- r*cos(eta[1])*sin(eta[2])
    theta[5]<- r*sin(eta[1])
    return(theta)
  }
  
  if(s==3) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])*cos(eta[3])
    theta[2]<- r*cos(eta[1])*cos(eta[2])*sin(eta[3])
    theta[3]<- r*cos(eta[1])*sin(eta[2])
    theta[4]<- r*sin(eta[1])
    return(theta)
    
  }
  else if(s==2) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])
    theta[2]<- r*cos(eta[1])*sin(eta[2])
    theta[3]<- r*sin(eta[1])
    return(theta)
  }
  
  else if(s==1) {
    theta[1]<- r*cos(eta[1])
    theta[2]<- r*sin(eta[1])
    return(theta)
  }
  
}


# The following function computes polar coordinates from the cartesian 
# Input :  x is the cartesian coordinate
# Outputs:  1) r  : the radius of the polar coordinates
#           2) eta: the (p-1) by 1 vector of polar coordinates

cart2polar <- function(x) {
  p =length(x)                  # number of covariates
  eta =  vector(length = p-1)   # vector for storing the radian angles
  r <- sqrt(sum(x^2))         # find the norm
  if(p==2) {
    eta[1]= atan(x[2]/x[1])
  }
  
  if(p==4) {
    eta[1] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[2] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[3] = atan(x[2]/x[1])
  }
  
  if(p==5) {
    eta[1] = atan(x[5]/sqrt(sum((c(x[1],x[2],x[3],x[4]))^2)))
    eta[2] = atan(x[4]/sqrt(sum((c(x[1],x[2],x[3]))^2)))
    eta[3] = atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[4] = atan(x[2]/x[1])
  }
  return(c(r,eta))
}

# Cost function W_n 
WnCost <- function(eta, X, tt, Y, kern, h) {
  theta <- polar2cart(eta, 1)      # cartesian coordinates from polar coordinates
  z <- X%*%theta        # creating the single index
  Yhat <- LocWassRegAMP(xin = z, qin = Y, optns= list(bwReg = h, qSup = tt, lower=20, 
                                                  upper=110, kernelReg= kern))  
  
  #get Wasserstein distance between the response and its prediction
  #pe <- matrix(0, nrow= nrow(Y), ncol= 1)
  #for (l in 1:nrow(Y)) {
  #  pe[l,] <- dist4den(list(x= tt, y= as.vector(Y_pred$qout[l,])), list(x= tt, y=as.numeric(Y[l,])), fctn_type = 'quantile')
  #}
  RSS <- sapply(1:nrow(Y), function(i) fdadensity:::trapzRcpp(X = tt, Y = (Y[i, ] - Yhat[i, ])^2))
  return(mean(RSS)) # compute the W_n function
}



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


