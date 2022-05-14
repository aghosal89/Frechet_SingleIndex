
## This is an alternative version of LocWassReg from the frechet package that 
## implements two changes:
##
## 1. The quadratic program is only set up if necessary.  So, if the local linear
##    fit is already increasing, no optimization will be attempted.
## 2. Includes an option to omit projection in case an extrinsic estimate is 
##    suitable for a given purpose. This can be specified via the option
##    'proj' in the optns list, which is set to TRUE by default.  Set to FALSE
##    to skip the projection step.


LocWassRegAMP <- function(xin, qin, xout = NULL, optns = list()) {
  
  if (!is.matrix(xin) & !is.vector(xin)) {
    stop("xin must be a matrix or vector")
  }
  if (is.vector(xin)) {
    xin <- matrix(xin, length(xin))
  }
  if (is.null(xout)) {
    xout <- xin
  }
  if (!is.matrix(xout) & !is.vector(xout)) {
    stop("xout must be a matrix or vector")
  }
  if (is.vector(xout)) {
    xout <- matrix(xout, length(xout))
  }
  if (ncol(xin) != ncol(xout)) {
    stop("xin and xout must have the same number of columns")
  }
  if (is.null(optns$bw)) {
    stop("optns$bw has no default values and must be input by user.")
  }
  if (!is.numeric(optns$bw) | (length(optns$bw) != ncol(xin))) {
    stop("optns$bw should be a numerical vector of length p.")
  }
  if (!is.matrix(qin) | !is.numeric(qin)) {
    stop("qin should be a numerical matrix.")
  }
  if (nrow(xin) != nrow(qin)) {
    stop("The number of rows of xin should be the same as the number of rows of qin.")
  }
  if (is.null(optns$ker)) {
    optns$ker <- "gauss"
  }
  if (is.null(optns$proj)) {
    optns$proj = TRUE
  }
  
  ker <- frechet:::kerFctn(optns$ker)
  K = function(x, h) {
    k = 1
    for (i in 1:p) {
      k = k * ker(x[, i]/h[i])
    }
    return(as.numeric(k))
  }
  k <- nrow(xout)
  n <- nrow(xin)
  m <- ncol(qin)
  p <- ncol(xin)
  getLFRweights = function(x0) {
    aux = K(xin - matrix(t(x0), nrow = n, ncol = length(x0), 
                         byrow = TRUE), optns$bw)
    mu0 = mean(aux)
    mu1 = colMeans(aux * (xin - matrix(t(x0), nrow = n, ncol = length(x0), 
                                       byrow = TRUE)))
    mu2 = 0
    for (i in 1:n) {
      mu2 = mu2 + aux[i] * (xin[i, ] - x0) %*% t(xin[i, 
      ] - x0)/n
    }
    sL = array(0, n)
    for (i in 1:n) {
      sL[i] = aux[i] * (1 - t(mu1) %*% solve(mu2) %*% (xin[i, 
      ] - x0))
    }
    s = sum(sL)
    return(sL/s)
  }
  A <- cbind(diag(m), rep(0, m)) + cbind(rep(0, m), -diag(m))
  if (!is.null(optns$upper) & !is.null(optns$lower)) {
    b0 <- c(optns$lower, rep(0, m - 1), -optns$upper)
  }
  else if(!is.null(optns$upper)) {
    A <- A[, -1]
    b0 <- c(rep(0, m - 1), -optns$upper)
  }
  else if (!is.null(optns$lower)) {
    A <- A[, -ncol(A)]
    b0 <- c(optns$lower, rep(0, m - 1))
  }
  else {
    A <- A[, -c(1, ncol(A))]
    b0 <- rep(0, m - 1)
  }
  Pmat <- as(diag(m), "sparseMatrix")
  Amat <- as(t(A), "sparseMatrix")
  qout <- sapply(1:k, function(j) {
    s = getLFRweights(xout[j,])
    s = as.vector(s)
    gx <- colMeans(qin * s) * n
    if(any(diff(gx) < 0) & optns$proj) {
      res <- do.call(osqp::solve_osqp, list(P = Pmat, q = -gx, 
                A = Amat, l = b0, pars = osqp::osqpSettings(verbose = FALSE)))
      gx <- sort(res$x)
    }
    return(gx)
  })
  qout = t(qout)
  return(qout)
}


