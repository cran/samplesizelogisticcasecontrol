#######################################
# Cont distribution X
#######################################
f3Cont <- function(x, alpha, beta, pdf){ 

  ret <- pdf(x)*px(x, alpha, beta)
  ret

} # END: f3Cont


fx1 <- function(x, alphastar, beta, prev, pdf) {
  pdf(x)*px(x,alphastar, beta)/prev #density of x in cases
}
fx0<-function(x, alphastar, beta, prev, pdf) {
  pdf(x)*(1-px(x,alphastar, beta))/(1-prev) #density x in controls
}
# Functions for computing the mean of x in the case-control population
f5 <- function(x, alpha, beta, pdf) {
  x*pdf(x)*px(x, alpha, beta) #x*f(x)*px
}
f6 <- function(x,alpha, beta, pdf) {
  x*pdf(x)*(1-px(x,alpha, beta)) #x*f(x)*(1-px)
}
 
# Function to compute the intercept in general population
getPopInt_dist_f <- function(beta, prev, pdf, pdf.range, interval=c(-7, 7), tol=0.0001) {

  a <- pdf.range[1]
  b <- pdf.range[2]

  pd <- function(alpha, beta, f3){
    #integrate(f3, a, b, alpha, beta, pdf)$value
    myintegrate_PDF(f3, a, b, method=2, alpha, beta, pdf)
  }
  f4 <- function(alpha, beta, prev, pd, f3){
    prev-pd(alpha, beta, f3) 
  }
  temp <- uniroot(f4, interval, tol=tol, beta, prev, pd, f3Cont)  
  root <- temp$root

  root

} # END: getPopInt_dist_f

# Function to compute the intercept in case-control population
getCCInt_dist_f <- function(alphastar, beta, prev, lam, pdf, pdf.range, f3Cont, 
                     interval=c(-7, 7), tol=0.0001) {

  a <- pdf.range[1]
  b <- pdf.range[2]

  fxcc <- function(x, alphastar, beta, lam, fx1, fx0, prev) {
    lam*fx1(x,alphastar,beta,prev, pdf)+(1-lam)*fx0(x,alphastar,beta,prev, pdf) #density of x in cc pop
  }
  fxccpx <- function(x, fxcc, alphastar, beta, lam, fx1, fx0, prev, px, alpha) {
    fxcc(x,alphastar, beta,lam,fx1,fx0,prev)*px(x,alpha,beta)  #fxcc*px
  }
  f7 <- function(alpha, fxccpx, fxcc, alphastar, beta, lam, fx1, fx0, prev, px) {
    #tmp <- integrate(fxccpx,a,b,fxcc,alphastar,beta,lam,fx1,fx0,prev,px,alpha)$value
    tmp <- myintegrate_PDF(fxccpx,a,b,method=2,fxcc,alphastar,beta,lam,fx1,fx0,prev,px,alpha)
    lam-tmp
  }

  temp <- uniroot(f7, interval, fxccpx, fxcc,alphastar,beta,lam,fx1,fx0,prev,px,tol=tol)
  root <- temp$root

  root

} # END: getCCInt_dist_f

# Compute information matrix
getInfo_dist_f <- function(alphacc, alphastar, beta, prev, lam, rangeF, pdf) {

  Intaa <- function(x,alphacc,alphastar,beta,prev,lam) {
    tempx <- px(x,alphacc,beta)
    temp  <- tempx*(1-tempx)
    ret   <- lam*fx1(x,alphastar,beta,prev,pdf)*temp +
             (1-lam)*fx0(x,alphastar,beta,prev,pdf)*temp
    ret
  }
  Intab <- function(x,alphacc,alphastar,beta,prev,lam) {
    tempx <- px(x,alphacc,beta)
    temp  <- tempx*(1-tempx)
    ret   <- lam*x*fx1(x,alphastar,beta,prev,pdf)*temp +
             (1-lam)*x*fx0(x,alphastar,beta,prev,pdf)*temp
    ret
  }
  Intbb <- function(x,alphacc,alphastar,beta,prev,lam) {
    tempx <- px(x,alphacc,beta)
    temp  <- tempx*(1-tempx)
    x2    <- x*x
    ret   <- lam*x2*fx1(x,alphastar,beta,prev,pdf)*temp +
             (1-lam)*x2*fx0(x,alphastar,beta,prev,pdf)*temp
    ret
  }
   
  a   <- rangeF[1]
  b   <- rangeF[2]
  Iaa <- myintegrate_PDF(Intaa, a,b, method=2,alphacc,alphastar,beta,prev,lam)
  Iab <- myintegrate_PDF(Intab, a,b, method=2,alphacc,alphastar,beta,prev,lam)
  Ibb <- myintegrate_PDF(Intbb, a,b, method=2,alphacc,alphastar,beta,prev,lam)

  info      <- matrix(0, 2, 2)
  info[1,1] <- Iaa
  info[1,2] <- Iab
  info[2,1] <- Iab
  info[2,2] <- Ibb

  info

} # END: getInfo_dist_f

# Variance
getVar_dist_f <- function(lam, info, Varx) {

  # NULL hyothesis
  c        <- lam*(1-lam)
  Var0beta <- 1/(c*Varx)
  Var0T    <- Varx/c

  # Alternative
  Ibba     <- info[2,2] - (info[1,2]*info[1,2])/info[1,1]
  Varabeta <- 1/Ibba
  VaraT    <- Ibba/c^2

  list(var.wald.null=Var0beta, var.wald.alt=Varabeta,
       var.score.null=Var0T, var.score.alt=VaraT)

} # END: getVar_dist_f

# Function to check the pdf
ss_check_pdf <- function(pdf, range) {

  x    <- seq(-10, 10, by=1)
  temp <- (x > range[1]) & (x < range[2])
  x    <- x[temp]
  y    <- pdf(x)
  temp <- (!is.finite(y)) | (y < 0)
  temp[is.na(temp)] <- TRUE
  n    <- sum(temp)
  if (n) {
    stop("ERROR: check distF and/or distF.support")
  }
 
  temp <- is.finite(range)
  if (any(temp)) {
    vec <- range[temp]
    for (v in vec) {
      temp <- pdf(v)
      if ((temp < 0) || (temp %in% c(NA, NaN))) {
        str <- paste("ERROR: ", v, " is not valid for distF.support", sep="")
        stop(str)
      }
    }
  }

  NULL

} # END: ss_check_pdf

######################################################################
# Main function for univariate case, continuous exposure, distribution
######################################################################
ss_uni_dist_f <- function(prev, logOR, pdf, pdf.range, size.2sided=0.05, 
                  power=0.9, lam=0.5, interval=c(-7, 7), tol=0.0001,
                  pdf.var=NULL, pdf.var.n=100000) {

  #  Check the pdf
  mypdf   <- NULL
  pdf.str <- check_pdf(pdf) 
  if (is.null(pdf.str)) pdf.str <- getObjName(pdf)
  if (length(pdf.range) != 2) stop("pdf.range must have length 2")

  # Create new function called mypdf and test it
  temp <- parsePDF(pdf.str, newPDF="mypdf", newRandF="myrpdf") 
  eval(parse(text=temp$pdf))
  ss_check_pdf(mypdf, pdf.range)

  # Save the str for generating random values 
  rpdf.str <- temp$rand
  varX     <- getVarX(pdf.var, rpdf.str, mypdf, n=pdf.var.n, lower=pdf.range[1],
                      upper=pdf.range[2]) 
  beta     <- logOR

  # Get the intercept in general population
  alphastar <- getPopInt_dist_f(beta, prev, mypdf, pdf.range, interval=interval, tol=tol) 

  # Function to compute the intercept in case-control population
  alphacc <- getCCInt_dist_f(alphastar, beta, prev, lam, mypdf, pdf.range, f3Cont, 
                     interval=interval, tol=tol)

  # Compute information matrix
  info <- getInfo_dist_f(alphacc, alphastar, beta, prev, lam, pdf.range, mypdf)

  a    <- pdf.range[1]
  b    <- pdf.range[2]
  #mux1 <- integrate(f5,a,b,alphastar,beta,mypdf)$value/prev  # mean x in cases
  #mux0 <- integrate(f6,a,b,alphastar,beta,mypdf)$value/(1-prev)  #mean x in controls
  mux1 <- myintegrate_PDF(f5,a,b,method=2,alphastar,beta,mypdf)/prev # mean x in cases
  mux0 <- myintegrate_PDF(f6,a,b,method=2,alphastar,beta,mypdf)/(1-prev) #mean x in controls

  # Get the variances for the Wald and score test
  temp           <- getVar_dist_f(lam, info, varX)
  var.wald.null  <- temp$var.wald.null
  var.wald.alt   <- temp$var.wald.alt
  var.score.null <- temp$var.score.null
  var.score.alt  <- temp$var.score.alt

  # Compute sample sizes for Wald test
  temp <- getSampleSize(power, size.2sided, logOR, var.wald.null, var.wald.alt)
  ss.wald.1 <- ceiling(temp$ss.eq1)
  ss.wald.2 <- ceiling(temp$ss.eq2)

  # Compute sample sizes for score test
  temp <- getSampleSize(power, size.2sided, mux1-mux0, var.score.null, var.score.alt)
  ss.score.1 <- ceiling(temp$ss.eq1)
  ss.score.2 <- ceiling(temp$ss.eq2)

  list(ss.wald.1=ss.wald.1, ss.wald.2=ss.wald.2, 
       ss.score.1=ss.score.1, ss.score.2=ss.score.2)

} # END: ss_uni_dist_f





##################
# Move to stat2.R
##################
# Function to check endpoint of integration
checkEndpoint_PDF <- function(a, fa) {

  a.vert    <- 0
  a.finite  <- is.finite(a) 
  fa.finite <- is.finite(fa)

  if (!is.numeric(a)) stop("Endpoint must be numeric")
  if (a.finite) {
    if ((fa.finite) && (fa < 0)) stop("PDF < 0 at endpoint")
    if (!fa.finite) a.vert <- 1
  } else {
    if (a %in% c(NA, NaN)) stop("Endpoint cannot be a missing value")
  }

  a.vert

} # END: checkEndpoint_PDF

getMid_PDF <- function(f, lower, upper, which, miny=10, maxy=100, ...) {

    v1 <- 10^(-(8:-4))
    if (which == 1) {
      xvec <- lower + v1
      temp <- xvec <= upper
      xvec <- xvec[temp]
    } else {
      xvec <- xvec - v1
      temp <- xvec >= lower
      xvec <- xvec[temp]
    }    
    yvec <- f(xvec, ...)

    temp <- (yvec >= miny) & (yvec <= maxy)
    n    <- sum(temp)
    if (n) {
      xvec <- xvec[temp]
      yvec <- yvec[temp]  
      if (which == 1) {
        # Lower
        ret <- xvec[n]
      } else {
        ret <- xvec[1]
      }
    } else {
      c   <- (miny + maxy)/2
      v   <- abs((yvec - c)*(yvec - c))
      i   <- which.min(v)
      ret <- xvec[i]
    }

    ret

  } # END: getMid


# Function to get midpoints of integration
getMidpoints_PDF <- function(f, lower, upper, ...) {

  # Check endpoints
  lower.finite <- is.finite(lower)
  upper.finite <- is.finite(upper)

  lower.vert   <- checkEndpoint_PDF(lower, f(lower, ...)) 
  upper.vert   <- checkEndpoint_PDF(upper, f(upper, ...)) 

  points <- lower
  if (lower.vert) {
    points <- c(points, getMid_PDF(f, lower, upper, 1, miny=10, maxy=100, ...))
  }
  if (upper.vert) {
    points <- c(points, getMid_PDF(f, lower, upper, 2, miny=10, maxy=100, ...))
  }
  points <- c(points, upper)

  if (length(points) == 2) {
    # Find "max" value
    v1   <- seq(-100, 100, by=0.01)
    temp <- (v1 > lower) & (v1 < upper)
    v1   <- v1[temp]
    if (!length(v1)) v1 <- NULL 
    v2   <- seq(101, 100000, by=1)
    temp <- (v2 > lower) & (v2 < upper)
    v2   <- v2[temp] 
    if (!length(v2)) v2 <- NULL 
    v3   <- seq(-100000, -101, by=1)
    temp <- (v3 > lower) & (v3 < upper)
    v3   <- v3[temp] 
    if (!length(v3)) v3 <- NULL 

    xvec   <- c(v3, v1, v2)
    yvec   <- f(xvec, ...)
    i      <- which.max(yvec)    
    points <- sort(c(points, xvec[i]))
  }
  
  points

} # END: getMidpoints_PDF

#  Function for general integration
myintegrate_PDF <- function(f, lower, upper, method=2, ...) {

  if (method == 1) {
    return(integrate(f, lower, upper, ...)$value)
  }

  # Get intervals of integration
  points  <- getMidpoints_PDF(f, lower, upper, ...)
  npoints <- length(points)
  nint    <- npoints - 1
  sum     <- 0
  for (i in 1:nint) {
    int <- integrate(f, points[i], points[i+1], ...)$value
    sum <- sum + int
  }

  sum

} # END: myintegrate_PDF

# Function to compute the mean of a cont rv X
getMean_PDF <- function(f, lower, upper, method=2, ...) {

  myfunc <- function(x, f2, ...) {
    x*f2(x, ...)
  }

  myintegrate_PDF(myfunc, lower, upper, method=method, f, ...)

} # END: getMean_PDF

# Function to compute the variance of a cont rv X
getVar_PDF <- function(f, lower, upper, method=2, EX=NULL, ...) {

  if (is.null(EX)) {
    EX <- getMean_PDF(f, lower, upper, method=method, ...)
  }

  myfunc <- function(x, f2, EX, ...) {
    (x - EX)*(x - EX)*f2(x, ...)
  }

  myintegrate_PDF(myfunc, lower, upper, method=method, f, EX, ...)

} # END: getVar_PDF




