# Function that uniroot function calls for computing intercept in general pop
#  fort the case of a ordinal exposure and no data
f1_dist_pmf <- function(alpha, beta, prev, p, xvals) {

  ret <- prev - sum(px(xvals, alpha, beta)*p)

  ret

} # END: f1_dist_pmf

getPopInt_dist_pmf <- function(beta, prev, p, xvals, interval=c(-7, 0), tol=0.0001) {

  temp <- uniroot(f1_dist_pmf, interval, tol=tol, beta, prev, p, xvals) 
  root <- temp$root

  root

} # END: getPopInt_dist_pmf

# Function to compute the prevalence of x=1 in the case-control population
pdcc_dist_pmf <- function(prev, alphastar, beta, p, xvals, lam) {

  tempx  <- px(xvals, alphastar, beta)
  #denom1 <- sum(p*tempx)
  #denom2 <- sum(p*(1-tempx))
  #px1    <- p*tempx/denom1
  #px0    <- p*(1-tempx)/denom2
  px1    <- p*tempx/prev
  px0    <- p*(1-tempx)/(1-prev)
  pxcc   <- lam*px1 + (1-lam)*px0

  list(px1.case=px1, px1.control=px0, px1=pxcc)

} # END: pdcc_dist_pmf

# Function uniroot function will call
f7_dist_pmf <- function(alpha, beta, lam, p, xvals){

  ret <- lam - sum(p*px(xvals, alpha, beta))
  ret

} # END: f7_dist_pmf

# Function to compute intercept in case-control population
getCCInt_dist_pmf <- function(beta, p, xvals, lam, interval=c(-7, 7), tol=0.0001) {

  temp <- uniroot(f7_dist_pmf, interval, tol=tol, beta, lam, p, xvals)
  root <- temp$root

  root

} # END: getCCInt_dist_pmf

# Information matrix
getInfo_dist_pmf <- function(alphacc, beta, lam, prev, p, xvals, px1, px0) {
  
  # p   Vector of Prob(X=i)
  # px1 Prevalence of X=i in cases
  # px0 Prevalence of X=i in controls   

  tempx <- px(xvals,alphacc,beta)
  val1  <- tempx*(1-tempx)
  pxv1  <- px1*val1
  pxv0  <- px0*val1
  Iaa   <- lam*sum(pxv1) + (1-lam)*sum(pxv0)
  Iab   <- lam*sum(pxv1*xvals) + (1-lam)*sum(pxv0*xvals)
  Ibb   <- lam*sum(pxv1*xvals*xvals) + (1-lam)*sum(pxv0*xvals*xvals)

  info      <- matrix(0, 2, 2)
  info[1,1] <- Iaa
  info[1,2] <- Iab
  info[2,1] <- Iab
  info[2,2] <- Ibb

  info

} # END: getInfo_dist_pmf

# Function to get the mean of a pmf
getMean_pmf <- function(xvals, px) {

  sum(xvals*px)

} # END: getMean_pmf

# Function to get the variance of a pmf
getVar_pmf <- function(xvals, px) {

  EX <- getMean_pmf(xvals, px)
  V  <- sum((xvals-EX)*(xvals-EX)*px)
  V

} # END: getVar_pmf

# Variance
getVar_dist_pmf <- function(px1, px0, xvals, lam, info) {
  
  c <- lam*(1-lam)

  #Varx <- getVar_pmf(xvals, p)
  #Var0beta <- 1/(c*Varx)
  #Var0T    <- Varx/c

  # NULL hyothesis
  temp     <- xvals*(lam*px1 + (1-lam)*px0)
  mu1bar   <- sum(temp)       # E(X)
  mu2bar   <- sum(xvals*temp) # E(X^2) 
  Iaa0     <- c
  Iab0     <- c*mu1bar
  Ibb0     <- c*mu2bar
  Ibba0    <- Ibb0-(Iab0*Iab0)/Iaa0
  Var0beta <- 1/Ibba0
  Var0T    <- Ibba0/c^2

  # Alternative
  Ibba     <- info[2,2] - (info[1,2]*info[1,2])/info[1,1]
  Varabeta <- 1/Ibba
  VaraT    <- Ibba/c^2

  list(var.wald.null=Var0beta, var.wald.alt=Varabeta,
       var.score.null=Var0T, var.score.alt=VaraT)

} # END: getVar_dist_pmf

###################################################################
# Main function for univariate case, ordinal exposure, distribution
###################################################################
ss_uni_dist_pmf <- function(prev, logOR, p, xvals, size.2sided=0.05, 
                     power=0.9, lam=0.5, interval=c(-7, 7), tol=0.0001) {

  beta <- logOR

  # Get the intercept in the general pop
  alphastar <- getPopInt_dist_pmf(beta, prev, p, xvals, interval=interval, tol=tol)

  # Compute the prevalence of x=1 in the case-control population
  temp <- pdcc_dist_pmf(prev, alphastar, beta, p, xvals, lam) 
  px1  <- temp$px1.case
  px0  <- temp$px1.control
  pxcc <- temp$px1

  # Get the intercept in the cc pop
  alphacc <- getCCInt_dist_pmf(beta, pxcc, xvals, lam, interval=interval, tol=tol)

  # Compute information matrix
  info <- getInfo_dist_pmf(alphacc, beta, lam, prev, p, xvals, px1, px0)

  # Get the variances for the Wald and score test
  temp           <- getVar_dist_pmf(px1, px0, xvals, lam, info)
  var.wald.null  <- temp$var.wald.null
  var.wald.alt   <- temp$var.wald.alt
  var.score.null <- temp$var.score.null
  var.score.alt  <- temp$var.score.alt

  # Compute sample sizes for Wald test
  temp <- getSampleSize(power, size.2sided, logOR, var.wald.null, var.wald.alt)
  ss.wald.1 <- ceiling(temp$ss.eq1)
  ss.wald.2 <- ceiling(temp$ss.eq2)

  # Compute sample sizes for score test
  mu1  <- getMean_pmf(xvals, px1) 
  mu0  <- getMean_pmf(xvals, px0) 
  temp <- getSampleSize(power, size.2sided, mu1-mu0, var.score.null, var.score.alt)
  ss.score.1 <- ceiling(temp$ss.eq1)
  ss.score.2 <- ceiling(temp$ss.eq2)

  list(ss.wald.1=ss.wald.1, ss.wald.2=ss.wald.2, 
       ss.score.1=ss.score.1, ss.score.2=ss.score.2)

} # END: ss_uni_dist_pmf

