# Function that uniroot function calls for computing intercept in general pop
#  for the case of a binary exposure and no data
f1_dist_bin <- function(alpha, beta, prev, p) {

  ret <- prev - p*px(1, alpha, beta) - (1-p)*px(0, alpha, beta)

  ret

} # END: f1_dist_bin

getPopInt_dist_bin <- function(beta, prev, p, interval=c(-7, 0), tol=0.0001) {

  temp <- uniroot(f1_dist_bin, interval, tol=tol, beta, prev, p) 
  root <- temp$root

  root

} # END: getPopInt_dist_bin

# Function to compute the prevalence of x=1 in the case-control population
pdcc_dist_bin <- function(prev, alphastar, beta, p, lam) {

  #temp0 <- px(0, alphastar, beta)
  temp1 <- px(1, alphastar, beta)
  #px1   <- p*temp1/(p*temp1 + (1-p)*temp0)  #in cases
  #px0   <- p*(1-temp1)/(p*(1-temp1) + (1-p)*(1-temp0))  #in controls
  px1   <- p*temp1/prev  #in cases
  px0   <- p*(1-temp1)/(1-prev)  #in controls
  pxcc  <- lam*px1 + (1-lam)*px0

  list(px1.case=px1, px1.control=px0, px1=pxcc)

} # END: pdcc_dist_bin

# Function uniroot function will call
f7_dist_bin <- function(alpha, beta, lam, p){

  ret <- lam - p*px(1, alpha, beta)-(1-p)*px(0, alpha, beta)
  ret

} # END: f7_dist_bin

# Function to compute intercept in case-control population
getCCInt_dist_bin <- function(beta, p, lam, interval=c(-7, 7), tol=0.0001) {

  temp <- uniroot(f7_dist_bin, interval, tol=tol, beta, lam, p)
  root <- temp$root

  root

} # END: getCCInt_dist_bin

# Function to compute intercept in case-control population
getCCInt_bin <- function(beta, obj, lam, prev, alphastar=NULL, 
                     interval=c(-7, 7), tol=0.0001) {

  # obj: matrix, probability or vector of probs
  if (is.matrix(obj)) {
    ret <- getCCInt_data(beta, obj, lam, prev, alphastar, interval=interval, tol=tol)
  } else {
    ret <- getCCInt_dist_bin(beta, obj, lam, interval=interval, tol=tol)
  }

  ret

} # END: getCCInt_bin

# Information matrix
getInfo_dist_bin <- function(alphacc, beta, lam, prev, p, px1, px0) {
  
  # p   Bernoulli parameter
  # px1 Prevalence of X=1 in cases
  # px0 Prevalence of X=1 in controls   

  temp1 <- px(1,alphacc,beta)
  temp0 <- px(0,alphacc,beta)
  Iaa   <- lam*(px1*temp1*(1-temp1)+(1-px1)*temp0*(1-temp0))+
           (1-lam)*(px0*temp1*(1-temp1)+(1-px0)*temp0*(1-temp0))
  Iab   <- lam*(px1*1*temp1*(1-temp1)+(1-px1)*0*temp0*(1-temp0))+
           (1-lam)*(px0*1*temp1*(1-temp1)+(1-px0)*0*temp0*(1-temp0))
  Ibb   <- lam*(px1*1*temp1*(1-temp1)+(1-px1)*0*temp0*(1-temp0))+
           (1-lam)*(px0*1*temp1*(1-temp1)+(1-px0)*0*temp0*(1-temp0))
 
  info      <- matrix(0, 2, 2)
  info[1,1] <- Iaa
  info[1,2] <- Iab
  info[2,1] <- Iab
  info[2,2] <- Ibb

  info

} # END: getInfo_dist_bin

# Variance
getVar_dist_bin <- function(px1, px0, lam, info) {
  
  c  <- lam*(1-lam)
  c2 <- c*c

  # NULL hyothesis
  mu1bar   <- lam*px1 + (1-lam)*px0
  mu2bar   <- mu1bar  #E(X)=E(X^2) because X=0 or 1
  Iaa0     <- c
  Iab0     <- c*mu1bar
  Ibb0     <- c*mu2bar
  Ibba0    <- Ibb0-(Iab0*Iab0)/Iaa0
  Var0beta <- 1/Ibba0
  Var0T    <- Ibba0/c2

  # Alternative
  Ibba     <- info[2,2] - (info[1,2]*info[1,2])/info[1,1]
  Varabeta <- 1/Ibba
  VaraT    <- Ibba/c2

  list(var.wald.null=Var0beta, var.wald.alt=Varabeta,
       var.score.null=Var0T, var.score.alt=VaraT)

} # END: getVar_dist_bin

###################################################################
# Main function for univariate case, binary exposure, distribution
###################################################################
ss_uni_dist_bin <- function(prev, logOR, p, size.2sided=0.05, power=0.9, lam=0.5,
                            interval=c(-7, 7), tol=0.0001) {

  beta <- logOR

  # Get the intercept in the general pop
  alphastar <- getPopInt_dist_bin(beta, prev, p, interval=interval, tol=tol)

  # Compute the prevalence of x=1 in the case-control population
  temp <- pdcc_dist_bin(prev, alphastar, beta, p, lam) 
  px1  <- temp$px1.case
  px0  <- temp$px1.control
  pxcc <- temp$px1

  # Get the intercept in the cc pop
  alphacc <- getCCInt_dist_bin(beta, pxcc, lam, interval=interval, tol=tol)

  # Compute information matrix
  info <- getInfo_dist_bin(alphacc, beta, lam, prev, p, px1, px0)

  # Get the variances for the Wald and score test
  temp           <- getVar_dist_bin(px1, px0, lam, info)
  var.wald.null  <- temp$var.wald.null
  var.wald.alt   <- temp$var.wald.alt
  var.score.null <- temp$var.score.null
  var.score.alt  <- temp$var.score.alt

  # Compute sample sizes for Wald test
  temp <- getSampleSize(power, size.2sided, logOR, var.wald.null, var.wald.alt)
  ss.wald.1 <- ceiling(temp$ss.eq1)
  ss.wald.2 <- ceiling(temp$ss.eq2)

  # Compute sample sizes for score test
  temp <- getSampleSize(power, size.2sided, px1-px0, var.score.null, var.score.alt)
  ss.score.1 <- ceiling(temp$ss.eq1)
  ss.score.2 <- ceiling(temp$ss.eq2)

  list(ss.wald.1=ss.wald.1, ss.wald.2=ss.wald.2, 
       ss.score.1=ss.score.1, ss.score.2=ss.score.2)

} # END: ss_uni_dist_bin

