


# Function to compute sample sizes from equations 1 and 2
getSampleSize <- function(gamma, size.2sided, theta, Var20, Var2a) {

  size    <- size.2sided/2  
  Z1gam   <- qnorm(gamma, mean=0, sd=1, lower.tail=TRUE, log.p=FALSE)
  Z1alpha <- qnorm(size, mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
  c1      <- (Z1gam+Z1alpha)^2

  # Equation 1
  ss1 <- ((Z1gam*sqrt(Var2a)+Z1alpha*sqrt(Var20))^2)/(theta)^2
  
  # Equation 2 
  ss2 <- c1*Var2a/(theta)^2  

  dim(ss1) <- NULL
  dim(ss2) <- NULL
  
  list(ss.eq1=ss1, ss.eq2=ss2)

} # END: getSampleSize 

# Function to compute D(D=1|x)
px <- function(x, alpha, beta)
{
  if (length(beta)==1) {
    u <- alpha + x*beta
  } else {
    u <- alpha + x%*%beta
  }

  ret <- 1/(1+exp(-u))
  ret

} # END: px

# Function to parse a pdf string
parsePDF <- function(pdf.name, newPDF="mypdf", newRandF="myrpdf") {

    pdf     <- removeWhiteSpace(pdf.name)
    vec     <- getVecFromStr(pdf, delimiter="(")
    fname   <- vec[1]
    len     <- nchar(pdf)
    # Remove last ")"
    if (substr(pdf, len, len) == ")") pdf <- substr(pdf, 1, len-1)
    vec     <- getVecFromStr(pdf, delimiter=",")
    temp    <- nchar(vec) > 0
    vec     <- vec[temp]
    vec     <- vec[-1]
    n       <- length(vec)
    pdf.str <- paste(fname, "(x", sep="")
    if (n) {
      temp    <- paste(vec, collapse=",", sep="")
      pdf.str <- paste(pdf.str, ",", temp, sep="") 
    }
    pdf.str <- paste(pdf.str, ")", sep="") 

    ret  <- paste(newPDF, " <- function(x) {", pdf.str, "}", sep="")

    pdf2 <- pdf.str
    substr(pdf2, 1, 1) <- "r"
    ret2 <- paste(newRandF, " <- function(x) {", pdf2, "}", sep="")

    list(pdf=ret, rand=ret2)

} # END: parsePDF

# Function to get the varX
getVarX <- function(varX, rpdf.str, pdf, n=100000, lower=-Inf, upper=Inf) {

  if (is.numeric(varX)) {
    if (length(varX) != 1) stop("ERROR with length(varX) > 1")
    if (varX <= 0) stop("ERROR with varX <= 0")
    return(varX)
  } 

  ret <- NULL
  if (is.null(varX)) {
    # Determine if there is a function available to generate random values
    if (n >= 10) {
      myrpdf <- NULL
      eval(parse(text=rpdf.str))
      ret <- try(var(myrpdf(n)), silent=TRUE)
      if ("try-error" %in% class(ret)) ret <- NULL  
    }  
    if (is.null(ret)) {
      # Integrate as a last resort
      ret <- getVar_PDF(pdf, lower, upper, method=2, EX=NULL)
    }    
  }

  ret

} # END: getVarX

# Function to check the pdf argument
check_pdf <- function(pdf) {

  if (is.function(pdf)) {
    # Test
    pdf.str <- NULL
  } else if (is.character(pdf)) {
    if (length(pdf) != 1) stop("pdf cannot have length > 1")
    pdf.str <- pdf
  } else {
    stop("pdf must be a function or a character string")
  }

  pdf.str

} # END: check_pdf

# Function top check for errors
checkForErrors <- function(prev, logOR, size.2sided, power, cc.ratio, interval,
                           tol, n.samples) {

  if ((prev <= 0) || (prev >= 1)) stop("ERROR: with prevalence")
  temp <- !is.finite(logOR)
  if (any(temp)) stop("ERROR: all logOR must be finite")
  if ((size.2sided <= 0) || (size.2sided >= 1)) stop("ERROR: size.2sided must be between 0 and 1")
  if ((power <= 0) || (power >= 1)) stop("ERROR: power must be between 0 and 1")
  if (length(interval) != 2) stop("ERROR: with interval")
  if ((tol <= 0) || (tol >= 1)) stop("ERROR: tol must be between 0 and 1")
  if ((cc.ratio <= 0) || (cc.ratio >= 1)) stop("ERROR: cc.ratio must be between 0 and 1")
  if (n.samples < 10) stop("ERROR: n.samples is too small")

  NULL

} # END: checkForErrors

# Function to get the data
getAndCheckData <- function(data, beta, min.nr=10) {

  nbeta    <- length(beta)
  dataFlag <- is.data.frame(data) || is.matrix(data)
  if (dataFlag) {
    temp <- rep(TRUE, nrow(data))
    for (i in 1:ncol(data)) {
      data[, i] <- as.numeric(data[, i])
      temp      <- temp & is.finite(data[, i])
    } 
    data <- data[temp, , drop=FALSE] 
  } else if (is.list(data)) {
    data.list <- check.file.list(data, op=NULL)
    data.list <- default.list(data.list, c("exposure"), list("ERROR"), error=c(1))
    data      <- loadData.table(data.list)
    covars    <- unique(data.list[["covars", exact=TRUE]])
    expvar    <- unique(data.list[["exposure", exact=TRUE]])
    ncovars   <- length(covars)
    nexpvar   <- length(expvar)
    if (nexpvar > 1) {
      # Interaction: 
      expvec <- rep(1, nrow(data))
      for (v in expvar) expvec <- expvec*as.numeric(data[, v])
    } else {
      expvec <- data[, expvar, drop=FALSE]
    }
    if (ncovars) {
      data <- data[, covars, drop=FALSE]
    } else {
      data <- NULL
    }
    data <- cbind(data, expvec)
    temp <- rep(TRUE, nrow(data))
    for (i in 1:ncol(data)) {
      data[, i] <- as.numeric(data[, i])
      temp      <- temp & is.finite(data[, i])
    } 
    data <- data[temp, , drop=FALSE] 
  } else{
    stop("ERROR: data must be a matrix, data frame or list")
  }

  data <- as.matrix(data)
  if (nrow(data) < min.nr) stop("ERROR: data contains too few rows")
  if (ncol(data) != nbeta) stop("ERROR: data does not contain correct number of columns")

  data

} # END: getAndCheckData

# Main function for binary exposure
sampleSize_binary <- function(prev, logOR, probXeq1=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, power=0.9, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001,
      n.samples=10000) {

  
  checkForErrors(prev, logOR, size.2sided, power, cc.ratio, interval,tol, n.samples)
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
    ret <- ss_uni_dist_bin(prev, logOR, probXeq1, size.2sided=size.2sided, power=power, 
                  lam=cc.ratio, interval=interval, tol=tol) 
  } else {
    # Covariates, data or data generator must be supplied
    if ((!distF.flag) && (!data.flag)) stop("ERROR: distF or data must be specified") 

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret <- ss_data(prev, logOR, data, size.2sided=size.2sided, power=power, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- ss_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                  power=power, lam=cc.ratio, interval=interval, tol=tol,
                  rpdf.n=n.samples)
    }
  }
  
  ret

} # END: sampleSize_binary

# Main function for ordinal exposure
sampleSize_ordinal <- function(prev, logOR, probX=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, power=0.9, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001,
      n.samples=10000) {

  
  checkForErrors(prev, logOR, size.2sided, power, cc.ratio, interval, tol, n.samples)
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
    ret <- ss_uni_dist_pmf(prev, logOR, probX, xvals, size.2sided=size.2sided, 
                     power=power, lam=cc.ratio, interval=interval, tol=tol)
  } else {
    # Covariates, data or data generator must be supplied
    if ((!distF.flag) && (!data.flag)) stop("ERROR: distF or data must be specified") 

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret <- ss_data(prev, logOR, data, size.2sided=size.2sided, power=power, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- ss_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                  power=power, lam=cc.ratio, interval=interval, tol=tol,
                  rpdf.n=n.samples)
    }
  }
  
  ret

} # END: sampleSize_ordinal

# Main function for continuous exposure
sampleSize_continuous <- function(prev, logOR, distF=NULL, distF.support=c(-Inf, Inf), 
      data=NULL, size.2sided=0.05, power=0.9, cc.ratio=0.5, interval=c(-100, 100), 
      tol=0.0001, n.samples=10000, distF.var=NULL) {

  checkForErrors(prev, logOR, size.2sided, power, cc.ratio, interval,tol, n.samples)
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

      ret <- ss_data(prev, logOR, data, size.2sided=size.2sided, power=power, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      ret <- ss_uni_dist_f(prev, logOR, distF, distF.support, size.2sided=size.2sided, 
                  power=power, lam=cc.ratio, interval=interval, tol=tol,
                  pdf.var=distF.var, pdf.var.n=n.samples)
    }
  } else {
    # Covariates, data or data generator must be supplied
    distF <- paste("rmvnorm(x, rep(0, ", nbeta, "))", sep="")

    if (data.flag) {
      # Get and check data
      data <- getAndCheckData(data, logOR, min.nr=10)

      ret  <- ss_data(prev, logOR, data, size.2sided=size.2sided, power=power, 
                     lam=cc.ratio, interval=interval, tol=tol) 
    } else {
      # Use the data generator
      ret <- ss_covars_dist_f(prev, logOR, rpdf=distF, size.2sided=size.2sided, 
                  power=power, lam=cc.ratio, interval=interval, tol=tol,
                  rpdf.n=n.samples)
    }
  }
  
  ret

} # END: sampleSize_continuous

# Function to compute sample size passing in data
sampleSize_data <- function(prev, logOR, data, size.2sided=0.05, power=0.9, 
                     cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) {

  checkForErrors(prev, logOR, size.2sided, power, cc.ratio, interval,tol, 1000)

  # Get and check data
  data <- getAndCheckData(data, logOR, min.nr=10)

  ret <- ss_data(prev, logOR, data, size.2sided=size.2sided, power=power, 
                     lam=cc.ratio, interval=interval, tol=tol) 
  ret

} # END: sampleSize_data

# Function to get intercept bounds
getIntBounds <- function(f, ...) {

  a  <- 10
  x0 <- -a
  x1 <- 0
  f0 <- f(x0, ...)
  f1 <- f(x1, ...)
  if (f0*f1 < 0) return(c(x0, x1))
  if (f1 == 0) return(c(-a, a))
  if (f0 == 0) return(c(-2*a, 0))

  x1 <- x1 + a
  for (i in 1:1000) {
    x0 <- x0*a
    x1 <- x1*a
    f0 <- f(x0, ...)
    f1 <- f(x1, ...)
    if (f0*f1 < 0) return(c(x0, x1))
  }
  stop("ERROR: Bounds for estimating the intercept could not be determined")


} # END: getIntBounds

