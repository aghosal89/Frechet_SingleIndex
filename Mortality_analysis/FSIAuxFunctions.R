
## This file contains auxiliary functions for use with FSIDenReg.R

# This is a function to compute the p-dimensional cartesian coordinates from 
# given polar coordinates in (p-1) dimension and radius.
# Inputs: 1) eta: a (p-1) dimensional vector
#         2) r : radius of the polar coordinates (set as default, r=1)

# codes below were written for p=2,3,4,5; which sufficed for our analysis.

polar2cart <- function(eta, r  = 1) {
  # read the length (p-1)
  s <- length(eta)
  # define vector for storing cartesian coordinates
  theta <- vector(length = s+1)
  
  # compute and return cartesian coordinates for a given p.
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
  if(s==2) {
    theta[1]<- r*cos(eta[1])*cos(eta[2])
    theta[2]<- r*cos(eta[1])*sin(eta[2])
    theta[3]<- r*sin(eta[1])
    return(theta)
  }
  
  if(s==1) {
    theta[1]<- r*cos(eta[1])
    theta[2]<- r*sin(eta[1])
    return(theta)
  }
}


# This function computes polar coordinates from the cartesian coordinates in
# p-dimensional euclidean space. 
# Input  :     x is the cartesian coordinate, p-dimensional vector of reals.
# Outputs:  1) r   : the radius of the polar coordinates.
#           2) eta : the vector of polar coordinates with length (p-1).

# This performs the reverse operation of the 'polar2cart' function above. This 
# function is written for p=2,3,4,5; which sufficed for our analysis.

cart2polar <- function(x) {
  # read dimension of the cartesian coordinate:
  p =length(x)      
  # define vector for storing polar coordinates of length (p-1)
  eta =  vector(length = p-1)   
  # computing the radius 
  r <- sqrt(sum(x^2))    
  
  # compute the polar coordinates for given p
  if(p==2) {
    eta[1]= atan(x[2]/x[1])
  }
  
  if(p==3) {
    eta[1]= atan(x[3]/sqrt(sum((c(x[1],x[2]))^2)))
    eta[2]= atan(x[2]/x[1])
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

# This function computes the Cost function for Fréchet Single Index model estimation
# Inputs: 1) eta  : polar coordinates.
#         2) X    : nxp matrix of predictors, n observations correspond to rows, p covariates to columns.
#         3) tt   : length m grid spanning [0, 1], used as grid for quantile functions.
#         4) Y    : nxm matrix of "observed" quantile functions on grid tt.
#         5) kern : string indicating the kernel to be used.
#         6) h    : bandwidth.

WnCost <- function(eta, X, tt, Y, kern, h) {
  # compute cartesian coordinates from polar coordinates
  theta <- polar2cart(eta, 1)      
  # compute the single index
  z <- X%*%theta       
  # compute the fréchet mean from the Wasserstein space
  Yhat <- LocWassRegAMP(xin = z, qin = Y, optns= list(bwReg = h, qSup = tt, lower=20, 
                                                  upper=110, kernelReg= kern))  
  
  #get array of Wasserstein distances between the response and its prediction
  RSS <- sapply(1:nrow(Y), function(i) fdadensity:::trapzRcpp(X = tt, Y = (Y[i, ] - Yhat[i, ])^2))
  # compute and return average distance as the cost function
  return(mean(RSS))
}


