
# function to produce the R_oplus_squared for Wasserstein space.
# need to read libraries 'frechet', 'fdadensity' to run this function

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
