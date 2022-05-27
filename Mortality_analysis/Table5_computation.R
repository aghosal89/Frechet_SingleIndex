

#######################################
# Computations for Table 5
# get Frechet R-squared for the models
#######################################

# set working directory~


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

source("frechet_Rsquared.R")

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
write.csv( df_mspe_folds,"MSPE_folds.csv")







########################
# As per referee inquiry, we run multiple regression models for the covariates

# Define the response
Y<- apply(quant_all, 1, mean)

fit1<- lm(Y~ X_ctr)
summary(fit1)$r.squared*100
sum(fit1$residuals^2)/(fit1$df.residual)

fit1.1<- lm(Y~ X_ctr[,1])
summary(fit1.1)$r.squared*100
sum(fit1.1$residuals^2)/(fit1.1$df.residual)

fit1.2<- lm(Y~ X_ctr[,2])
summary(fit1.2)$r.squared*100
sum(fit1.2$residuals^2)/(fit1.2$df.residual)

fit1.3<- lm(Y~ X_ctr[,3])
summary(fit1.3)$r.squared*100
sum(fit1.3$residuals^2)/(fit1.3$df.residual)

fit1.4<- lm(Y~ X_ctr[,4])
summary(fit1.4)$r.squared*100
sum(fit1.4$residuals^2)/(fit1.4$df.residual)

fit1.5<- lm(Y~ X_ctr[,5])
summary(fit1.5)$r.squared*100
sum(fit1.5$residuals^2)/(fit1.5$df.residual)


# Single Index model

# This function computes the Cost function for regression Single Index model estimation
# Inputs: 1) eta  : polar coordinates.
#         2) X    : nxp matrix of predictors, n observations correspond to rows, p covariates to columns.
#         4) Y    : nxm matrix of "observed" responses.

WnCost_2 <- function(eta, X, Y) {
  # compute cartesian coordinates from polar coordinates
  theta <- polar2cart(eta, 1)      
  # compute the single index
  z <- X%*%theta       
  # compute the fréchet mean from the Wasserstein space
  fit_ <- lm(Y~z)
  
  mse <- (sum(fit_$residuals^2))/fit_$df.residual
  # compute and return average distance as the cost function
  return(mse)
}




## This function fits the SI index model in the case of responses in real line.
## See file FSIAUXFunctions.R for auxiliary functions used within this script.  
## 
##
## Inputs:
##
## X       - nxp matrix of predictors, n observations correspond to rows, p variables to columns
## Y       - nxm matrix of "observed" quantile functions on grid tt
## nsp     - integer giving the number of starting points in each dimension to be 
##           used by optim. A lattice of points will be created by constructing 
##           an equally spaced grid for each of the (p - 1) hyperspherical coordinates
##           used to represent theta in the optimization. Default is 3
## L       - an integer specifying number of starting points to use.  If L = 0 (default),
##           all of the starting points in the lattice will be utilized.  Otherwise,
##           L of these will be randomly and uniformly sampled
##
## Output: List with the following elements
##
## thetaHat - length p vector giving the estimated coefficient
## fnvalue  - achieved minimum value of the criterion function for estimating theta
## etaStart - matrix with (p - 1) columns, each row indicating a unique starting value
##            used in optimization for estimating theta
## optInf   - list containing information about optimization routine for each
##            starting point



SIReg <- function(X = NULL, Y = NULL, nsp = 3){
  
  # Perform checks
  if(is.null(X) | !is.matrix(X)){
    stop('Must provide covariates as a matrix X')
  }
  if(is.null(Y) | !is.matrix(Y)){
    stop('Must provide quantile functions as a matrix Y')
  }
  if(nrow(Y) != nrow(X)){
    stop('Dimensions of Y do not match X inputs')
  }
  if(!is.numeric(nsp) | length(nsp) != 1 | mod(nsp, 1) != 0 | nsp <= 0){
    message('Invalid specification of input nsp, resetting to default')
    nsp <- 3
  }
  
  # Create grid of starting values for optimization
  
  # compute dimension of covariate
  p <- ncol(X)
  # specified spacing between staring points in each coordinate
  spc <- pi/nsp  
  # equally spaced starting points in polar coordinates
  f <- lapply(1:(p - 1), function(j) seq(-pi/2 + spc/2, pi/2 - spc/2, by = spc)) 
  # create grid of starting values
  etaStart <- as.matrix(expand.grid(f))   
  
  ## To provide information about optimization as output
  optInf <- list()
  
  # provide criteria for termination of the algorithm
  optim_optns <- list(factr = 1e11, maxit = 100)
  
  WnMin <- rep(NA, nrow(etaStart))
  etaMin <- matrix(NA, nrow = nrow(etaStart), ncol = p - 1)
  
  # main optimization loop over starting values
  for (k in 1:nrow(etaStart)) {
    WnOpt <- optim(par = etaStart[k, ], fn = WnCost_2, method = "L-BFGS-B",
                   lower = -pi/2, upper = pi/2, control = optim_optns,
                   X = X, Y = Y)
    
    optInf[[k]] <- WnOpt
    
    WnMin[k] <- WnOpt$value
    etaMin[k, ] <- WnOpt$par
  }
  
  # the optimizer, i.e. thetaHat    
  thetaHat <- polar2cart(etaMin[which.min(WnMin),], 1)
  
  optvalue <- min(WnMin)  # the minimized cost
  
  fit_ <- lm(Y ~ X %*% thetaHat)
  Yout <- as.vector(fit_$fitted.values)
  # Create return variable as a list
  return(list(thetaHat = thetaHat, fnvalue=optvalue, Yout = Yout, 
              etaStart <- etaStart, optInf = optInf))
  
}

library("numbers")
Y<- as.matrix((apply(quant_all, 1, mean)), n,1)
sd=SIReg(X = X_ctr, Y = Y, nsp = 3)


sd$thetaHat
sd$fnvalue

