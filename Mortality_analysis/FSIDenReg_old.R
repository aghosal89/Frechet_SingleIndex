
## This function fits the FSI index model in the case of distributional responses 
## for a given bandwidth and kernel function.  See file FISAUXFunctions.R for 
## auxiliary functions used within this script.  This file must be sourced prior
## to using FSIDenReg.  The frechet and pracma packages must also be loaded
##
## Inputs:
##
## X       - nxp matrix of predictors, observations correspond to rows, variables to columns
## tt      - length m grid spanning [0, 1], used as grid for quantile functions
## Y       - nxm matrix of "observed" quantile functions on grid tt
## h       - bandwidth for smoothing
## kern    - string indicating the kernel to be used (options are 'gauss' or 'epan'),
##           default is 'gauss'
## Xout    - qxp matrix of predictors at which fitted quantile functions will be 
##           produced, defaults to NA
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
## Yout     - if Xout is provided, a qxm matrix of fitted quantile functions; 
##            otherwise, it will be empty
## optInf   - list containing information about optimization routine for each
##            starting point


# Note, to compute the value W_n(theta) for a given theta, use
# locDenReg(xin = X%*%theta, qin = Y, xout = Xout%*%theta, 
#           options = list(bwReg = h, kernelReg = kern, qSup = tt))
# Then compute squared Wasserstein distances between Y and the fitted
# quantile functions in qout to get W_n for that theta



FSIDenReg <- function(X = NULL, tt = NULL, Y = NULL, h = NULL, kern = 'gauss', 
                      Xout = NULL, nsp = 3, L = 0){
  
    # Perform checks
    if(is.null(X) | !is.matrix(X)){
      stop('Must provide covariates as a matrix X')
    }
    if(is.null(tt)){
      stop('Must provide grid vector tt for quantil functions')
    }
    if(is.null(Y) | !is.matrix(Y)){
      stop('Must provide quantile functions as a matrix Y')
    }
    if(nrow(Y) != nrow(X) | ncol(Y) != length(tt)){
      stop('Dimensions of Y do not match tt and/or X inputs')
    }
    if(!is.numeric(h) | length(h) != 1 | h <= 0){
      stop('Bandwidth h must be a positive scalar')
    }
    if(!(kern %in% c('rect', 'gauss', 'epan', 'gausvar', 'quar'))){
      stop('Invalid kernel choice')
    }
    if(!is.null(Xout) & (!is.matrix(Xout) | ncol(Xout) != ncol(X))){
      stop('Xout must be a matrix with number of columns matching that of X')
    }
    if(!is.numeric(nsp) | length(nsp) != 1 | mod(nsp, 1) != 0 | nsp <= 0){
      message('Invalid specification of input nsp, resetting to default')
      nsp <- 3
    }
    if(!is.numeric(L) | length(L) != 1 | mod(L, 1) != 0 | L < 0){
      message('Invalid specification of input L, resetting to default')
      nsp <- L
    }
  
    # Create grid of starting values for optimization

    p <- ncol(X)   # dimension of covariate
    spc <- pi/nsp  # specified spacing between staring points in each coordinate
    f <- seq(-pi/2 + spc/2, pi/2 - spc/2, by = spc)  # equally spaced starting points in polar coordinates
    etaStart <- as.matrix(expand.grid(f,f,f,f))  # create grid of starting values 
    if(L != 0) {
      smp <- sample.int(nrow(etaStart), min(L, nrow(etaStart)))
      etaStart <- etaStart[smp,]
    }
    
    #### For exploration only, this will be deleted
    optInf <- list()
    ####
    
    optim_optns = list(factr = 1e11, maxit = 100)
    #fopt <- matrix(0, nrow = nrow(etaStart), ncol= 1)
    #xopt<- matrix(0, nrow = nrow(etaStart), ncol= p - 1)
    WnMin <- rep(NA, nrow(etaStart))
    etaMin <- matrix(NA, nrow = nrow(etaStart), ncol = p - 1)
    
    
    for (k in 1:nrow(etaStart)) {
      WnOpt <- optim(par = etaStart[k, ], fn = WnCost, method = "L-BFGS-B",
                     lower = -pi/2, upper = pi/2, control = optim_optns,
                     X = X, tt = tt, Y = Y, kern = kern, h = h)
      
      optInf[[k]] <- WnOpt
      #opt <- optimset(TolFun= 1e-4) # decide on tolerance for the algorithm
      #run <- fminbnd(nm_optfun, as.vector(start[k,]), 
      #               c(-pi/2,-pi/2,-pi/2,-pi/2), 
      #               c(pi/2,pi/2,pi/2,pi/2), opt)
      #fopt[k]<- as.vector(neldermead.get(this=run, key ='fopt'))
      #xopt[k,]<- as.vector(neldermead.get(this=run, key ='xopt'))
      WnMin[k] <- WnOpt$value
      etaMin[k, ] <- WnOpt$par
    }
    
    # the optimizer, i.e. thetaHat    
    #theta_hat<- polar2cart(xopt[which.min(fopt),],1)
    thetaHat <- polar2cart(etaMin[which.min(WnMin),], 1)
    
    optvalue <- min(WnMin)
    # Lastly, if Xout is not empty, compute predictions using LocDenReg
    Yout <- NA
    if(!is.null(Xout)){
      Yout <- LocWassRegAMP(xin = X%*%thetaHat, qin = as.matrix(Y), 
                     xout = Xout%*%thetaHat, 
                  optns = list(bwReg = h, kernelReg= kern, lower = 20, 
                               upper = 110, qSup = tt))
    }
    # Create return variable as a list
    return(list(thetaHat = thetaHat, fnvalue=optvalue, Yout = Yout, optInf = optInf))
    
}


