

# Function to compute power from equations 1 and 2
getPower <- function(sampleSize, size.2sided, theta, Var20, Var2a) {

  size    <- size.2sided/2 
  Z1alpha <- qnorm(size, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)  
  c1      <- sqrt(abs(sampleSize*theta*theta/Var2a))

  # Equation 1
  Zgam1  <- c1 - Z1alpha*sqrt(Var20/Var2a)
  pow1   <- pnorm(Zgam1, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)  

  # Equation 2 
  Zgam2 <- c1 - Z1alpha  
  pow2  <- pnorm(Zgam2, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)  

  dim(pow1) <- NULL
  dim(pow2) <- NULL
  
  list(pow.eq1=pow1, pow.eq2=pow2)

} # END: getPower

# Main function for binary exposure
power_binary <- function(prev, logOR, probXeq1=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, sampleSize=1000, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) {

  n.samples <- sampleSize  
  checkForErrors(prev, logOR, size.2sided, NULL, cc.ratio, interval,tol, n.samples, sampleSize=sampleSize)
  data.flag  <- !is.null(data)
  p.flag     <- !is.null(probXeq1)
  if (p.flag) {
    if ((probXeq1 <= 0) || (probXeq1 >= 1)) stop("ERROR: probXeq1 must be between 0 and 1")
  }
  distF.flag <- !is.null(distF)
  
  
  nbeta <- length(logOR)
  if (nbeta == 1) {
    # No covariates
    if (!p.flag) {
      if (!data.flag) stop("ERROR: probXeq1 or data must be specified")  
      
      # Get and check the data
      data     <- getAndCheckData(data, logOR, min.nr=10)
      vec      <- data[, ncol(data)]
      if (!all(vec %in% 0:1)) stop("ERROR: binary exposure must be codeds as 0-1")

      probXeq1 <- sum(vec %in% 1)/length(vec)
    }       
    ret <- pow_uni_dist_bin(prev, logOR, probXeq1, size.2sided=size.2sided, sampleSize=sampleSize, 
                  lam=cc.ratio, interval=interval, tol=tol) 
  } else {
    # Covariates, data or data generator must be supplied
    if ((!distF.flag) && (!data.flag)) stop("ERROR: distF or data must be specified") 

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret <- pow_data(prev, logOR, data, size.2sided=size.2sided, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- pow_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                  sampleSize=sampleSize, lam=cc.ratio, interval=interval, tol=tol)
    }
  }
  
  ret

} # END: power_binary

# Main function for ordinal exposure
power_ordinal <- function(prev, logOR, probX=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, sampleSize=1000, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) {

  n.samples <- sampleSize
  checkForErrors(prev, logOR, size.2sided, NULL, cc.ratio, interval, tol, n.samples, sampleSize=sampleSize)
  data.flag  <- !is.null(data)
  p.flag     <- !is.null(probX)
  if (p.flag) {
    temp <- (probX <= 0) | (probX >= 1)
    temp[!is.finite(temp)] <- TRUE
    if (any(temp)) stop("ERROR: probX must be between 0 and 1")
    if (sum(probX) != 1) stop("ERROR: probX must sum to 1")
  }
  distF.flag <- !is.null(distF)
  
  nbeta <- length(logOR)
  if (nbeta == 1) {
    # No covariates
    if (!p.flag) {
      if (!data.flag) stop("ERROR: probX or data must be specified")  
      
      # Get and check the data
      data     <- getAndCheckData(data, logOR, min.nr=10)

      vec      <- data[, ncol(data)]
      n        <- length(vec)
      uvec     <- sort(unique(vec))
      nvec     <- length(uvec)
      probX    <- rep(NA, nvec)
      for (i in 1:nvec) {
        probX[i] <- sum(vec %in% uvec[i])/n
      }
      xvals <- 0:(nvec-1)
      if (!all(uvec == xvals)) {
        cat("\nNOTE: ordinal exposure variable X in data did not have values 0, 1, ..., k\n")
        cat("X is being re-coded as:\n")
        temp <- cbind(uvec, xvals)
        colnames(temp) <- c("X.original", "X")
        print(temp)  
      } 
    } else {
      xvals <- 0:(length(probX)-1)
    }       
    ret <- pow_uni_dist_pmf(prev, logOR, probX, xvals, size.2sided=size.2sided, 
                    sampleSize=sampleSize, lam=cc.ratio, interval=interval, tol=tol)
  } else {
    # Covariates, data or data generator must be supplied
    if ((!distF.flag) && (!data.flag)) stop("ERROR: distF or data must be specified") 

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret <- pow_data(prev, logOR, data, size.2sided=size.2sided, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- pow_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                 sampleSize=sampleSize, lam=cc.ratio, interval=interval, tol=tol)
    }
  }
  
  ret

} # END: power_ordinal

# Main function for continuous exposure
power_continuous <- function(prev, logOR, distF=NULL, distF.support=c(-Inf, Inf), 
      data=NULL, size.2sided=0.05, sampleSize=1000, cc.ratio=0.5, interval=c(-100, 100), 
      tol=0.0001, distF.var=NULL) {

  n.samples <- sampleSize
  checkForErrors(prev, logOR, size.2sided, NULL, cc.ratio, interval,tol, n.samples, sampleSize=sampleSize)
  data.flag  <- !is.null(data)
  distF.flag <- !is.null(distF)
  nbeta      <- length(logOR)

  if (nbeta == 1) {
    if (!distF.flag) {
      # Default to N(0,1)
      distF     <- dnorm
      distF.var <- 1
    }
    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret <- pow_data(prev, logOR, data, size.2sided=size.2sided, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      ret <- pow_uni_dist_f(prev, logOR, distF, distF.support, size.2sided=size.2sided, 
                  sampleSize=sampleSize, lam=cc.ratio, interval=interval, tol=tol,
                  pdf.var=distF.var)
    }
  } else {
    # Covariates, data or data generator must be supplied
    distF <- paste("rmvnorm(x, rep(0, ", nbeta, "))", sep="")

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret  <- pow_data(prev, logOR, data, size.2sided=size.2sided, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- pow_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                  sampleSize=sampleSize, lam=cc.ratio, interval=interval, tol=tol)
    }
  }
  
  ret

} # END: power_continuous

# Function to compute power passing in data
power_data <- function(prev, logOR, data, size.2sided=0.05, sampleSize=1000,
                     cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) {

  checkForErrors(prev, logOR, size.2sided, NULL, cc.ratio, interval,tol, 1000, sampleSize=sampleSize)

  # Get and check data
  data <- getAndCheckData(data, logOR, min.nr=10)

  ret <- pow_data(prev, logOR, data, size.2sided=size.2sided, sampleSize=sampleSize,
                     lam=cc.ratio, interval=interval, tol=tol) 
  ret

} # END: power_data
