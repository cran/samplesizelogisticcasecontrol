###############################################################################
# To average with respect to the distribuion in the general population, 
#   average over x in controls (or sample from gen pop)
# To average with respect to the distribution in cases, 
#   average over x with weight px(x)/prev in controls (or sample from gen pop)
# To average with respect to the distribution in controls, average over x 
#   with weight (1-px(x))/(1-prev) in controls (or sample from gene pop)
###############################################################################

# Function to compute average px using pilot data
avepx_data <- function(alpha, beta, data) {
  # Data: matrix of data corresponding to beta 

  pxvec <- px(data, alpha, beta)
  avepx <- mean(pxvec, na.rm=TRUE)

  avepx

} # END: avepx_data

# Function that uniroot function calls for computing intercept in general pop
f1_data <- function(alpha, beta, prev, data) {

  ret <- prev - avepx_data(alpha, beta, data)

  ret

} # END: f1_data

# Function to compute the intercept in the general population with 
#  disease prevalence prev using pilot data 
getPopInt_data <- function(beta, prev, data, interval=c(-7, 0), tol=0.0001) {

  temp <- uniroot(f1_data, interval, tol=tol, beta, prev, data) 
  root <- temp$root

  root

} # END: getPopInt_data

#####################################################################
#  Now compute intercept alpha in case-control population
#####################################################################
pdcc_data <- function(alpha, beta, data, lam, prev, alphastar) {

  pxvec   <- px(data, alphastar,beta) #prob disease in gen pop
  pxccvec <- px(data, alpha, beta)  #prob disease in cc population
  px2     <- pxccvec*pxvec  #for average over cases
  px3     <- (1-pxvec)*pxccvec #for average over controls
  pdccvec <- lam*px2/prev+(1-lam)*px3/(1-prev)
  pdcc    <- mean(pdccvec, na.rm=TRUE)

  pdcc

} # END: pdcc_data

# Function uniroot function will call
f7_data <- function(alpha, beta, data, lam, prev, alphastar) {

  ret <- lam - pdcc_data(alpha, beta, data,lam, prev, alphastar)
  ret

} # END: f7_data

# Function to compute intercept in case-control population
getCCInt_data <- function(beta, data, lam, prev, alphastar, 
                     interval=c(-7, 7), tol=0.0001) {

  temp <-uniroot(f7_data, interval, beta, data, lam, prev, alphastar, tol=tol)
  root <- temp$root

  root

} # END: getSourceInt_data

####################################################################################
# For score tests
#####################################################################################
# Function to compute the null scores
nullscore <- function(ab0, beta, data, lam, prev, alphastar) {

  # Last column of data will be for variable of interest
  
  nc        <- ncol(data)
  ncov      <- nc - 1
  nullscore <- rep(0, nc)
  pxgen     <- px(data,alphastar,beta) #prob disease in gen pop
  if (ncov) {
    x1   <- data[, 1:ncov, drop=FALSE]
    pxcc <- px(x1, ab0[1], ab0[-1])  #prob disease in cc population
  } else {
    pxcc <- ab0[1]
  }

  temp1        <- lam*(1-pxcc)*pxgen/prev
  temp2        <- (1-lam)*pxcc*(1-pxgen)/(1-prev)
  u1           <- temp1 - temp2
  nullscore[1] <- mean(u1, na.rm=TRUE)

  if (ncov) {
    for (i in 1:ncov) {
      u2 <- lam*x1[, i]*(1-pxcc)*pxgen/prev - (1-lam)*x1[, i]*pxcc*(1-pxgen)/(1-prev)
      nullscore[i+1] <- mean(u2, na.rm=TRUE)
    }
  }

  nullscore

} # END: nullscore

# Functions to solve for the cc intercept and beta1:betam under the null: beta_m+1=0
scorenorm <- function(ab0, beta, data, lam, prev, alphastar) {

  nullscores <- nullscore(ab0, beta, data, lam, prev, alphastar)
  scorenorm  <- sum(nullscores*nullscores)

  scorenorm

} # END: scorenorm

# Compute null estimates of intercept and slope in cc sample
nlm_data <- function(f, initVal, beta, data, lam, prev, alphastar) {

  ab0 <- nlm(f,initVal,beta,data,lam,prev,alphastar)$estimate
  ab0  

} # END: nlm_data

# Compute the information matrix under the full (alternative) model
getInfo_data <- function(alpha,beta,data,lam,prev,alphastar) {
  
  nbeta <- length(beta)
  nv    <- nbeta + 1
  Info  <- matrix(rep(0, nv*nv), nrow=nv, ncol=nv)
  pxgen <- px(data,alphastar,beta) #prob disease in gen pop
  pxcc  <- px(data,alpha,beta)  #prob disease in cc population

  tmp1 <- lam*pxcc*(1-pxcc)*pxgen/prev
  tmp2 <- (1-lam)*pxcc*(1-pxcc)*(1-pxgen)/(1-prev) 
  tmp  <- tmp1 + tmp2

  # For intercept
  Info[1,1] <- mean(tmp, na.rm=TRUE)

  # For intercept and other vars
  for (i in 1:nbeta) {
    col         <- i + 1
    temp        <- data[, i]
    vec         <- temp*tmp
    Info[1,col] <- mean(vec, na.rm=TRUE)
    Info[col,1] <- Info[1, col]
  }

  # For other vars
  for (i in 1:nbeta) {
    row  <- i + 1
    v1   <- data[, i]
    for (j in 1:nbeta) {
      col           <- j + 1
      v2            <- data[, j]
      vec           <- v1*v2*tmp
      Info[row,col] <- mean(vec, na.rm=TRUE)
      Info[col,row] <- Info[row, col]
    }
  }

  Info

} # END: getInfo_data

# Compute the information matrix at the null
getInfo0_data <- function(ab0,beta,data,lam,prev,alphastar) {
  
  nbeta  <- length(beta)
  nv     <- nbeta + 1
  ncov   <- nbeta - 1
  Info   <- matrix(rep(0, nv*nv), nrow=nv, ncol=nv)
  pxgen  <- px(data,alphastar,beta) #prob disease in gen pop
  alpha0 <- ab0[1]
  nc     <- ncol(data)
  ncm1   <- nc - 1
  if (ncov) {
    beta0 <- ab0[-1]
    pxcc  <- px(data[, 1:ncm1, drop=FALSE],alpha0,beta0)  #prob disease in cc population
  } else {
    pxcc  <- alpha0
  }
  
  tmp1 <- lam*pxcc*(1-pxcc)*pxgen/prev
  tmp2 <- (1-lam)*pxcc*(1-pxcc)*(1-pxgen)/(1-prev) 
  tmp  <- tmp1 + tmp2

  # For intercept
  Info[1,1] <- mean(tmp, na.rm=TRUE)

  # For intercept and other vars
  for (i in 1:nbeta) {
    col         <- i + 1
    temp        <- data[, i]
    vec         <- temp*tmp
    Info[1,col] <- mean(vec, na.rm=TRUE)
    Info[col,1] <- Info[1, col]
  }

  # For other vars
  for (i in 1:nbeta) {
    row  <- i + 1
    v1   <- data[, i]
    for (j in 1:nbeta) {
      col           <- j + 1
      v2            <- data[, j]
      vec           <- v1*v2*tmp
      Info[row,col] <- mean(vec, na.rm=TRUE)
      Info[col,row] <- Info[row, col]
    }
  }

  Info

} # END: getInfo_data

# Function to compute the expected value of the score for beta2 at null 
#   values ab0=(alpha0,beta10) and the null variance
U2score_data <- function(ab0,beta,data,lam,prev,alphastar) {

  nc    <- ncol(data)
  ncm1  <- nc - 1
  pxgen <- px(data,alphastar,beta) #prob disease in gen pop
  if (ncm1) {
    pxcc <- px(data[, 1:ncm1, drop=FALSE],ab0[1],ab0[-1])  #prob disease in cc population
  } else {
    pxcc <- ab0[1]
  }
  x2      <- data[, nc]
  U2s     <- lam*x2*(1-pxcc)*pxgen/prev - (1-lam)*x2*pxcc*(1-pxgen)/(1-prev)
  U2score <- mean(U2s, na.rm=TRUE)

  U2score

} # END: U2score_data

# Function for computing the variance
getVar_data <- function(info0, infoa) {
  
  nv        <- ncol(info0)
  nvm1      <- nv - 1 
  I11       <- info0[1:nvm1, 1:nvm1, drop=FALSE]
  I12       <- info0[1:nvm1, nv, drop=FALSE]
  I22       <- info0[nv, nv]
  VarU20    <- I22 - t(I12)%*%solve(I11)%*%I12
  Varbeta20 <- 1/VarU20

  I11       <- infoa[1:nvm1, 1:nvm1, drop=FALSE]
  I12       <- infoa[1:nvm1, nv, drop=FALSE]
  I22       <- infoa[nv, nv]
  VarU2a    <- I22 - t(I12)%*%solve(I11)%*%I12
  Varbeta2a <- 1/VarU2a

  list(var.wald.null=Varbeta20, var.wald.alt=Varbeta2a,
       var.score.null=VarU20, var.score.alt=VarU2a)

} # END: getVar_data

###################################################################
# Main function for pilot data, covariates
###################################################################
ss_data <- function(prev, beta, data, size.2sided=0.05, power=0.9, lam=0.5,
                    interval=c(-7, 7), tol=0.0001) {

  nbeta <- length(beta)
  ncov  <- nbeta - 1
  data  <- as.matrix(data)
 
  # Compute the intercept in the general population 
  alphastar <- getPopInt_data(beta, prev, data, interval=interval, tol=tol)

  # Compute intercept alpha in case-control population
  alphacc <- getCCInt_data(beta, data, lam, prev, alphastar, 
                     interval=interval, tol=tol)

  # Compute null estimates of intercept and slope in cc sample
  initVal <- rep(0, nbeta)
  ab0     <- nlm_data(scorenorm, initVal, beta, data, lam, prev, alphastar)

  # Compute information matrix under the alternative
  infoa <- getInfo_data(alphacc,beta,data,lam,prev,alphastar)

  # Compute information matrix under the null
  info0 <- getInfo0_data(ab0,beta,data,lam,prev,alphastar)

  # Compute the mean, theta, of the U2 score for beta2 under the alternative
  U20 <- U2score_data(ab0,beta,data,lam,prev,alphastar)

  # Get the variances for the Wald and score test
  temp           <- getVar_data(info0, infoa)

  var.wald.null  <- temp$var.wald.null
  var.wald.alt   <- temp$var.wald.alt
  var.score.null <- temp$var.score.null
  var.score.alt  <- temp$var.score.alt

  # Compute sample sizes for Wald test
  temp <- getSampleSize(power, size.2sided, beta[nbeta], var.wald.null, var.wald.alt)
  ss.wald.1 <- ceiling(temp$ss.eq1)
  ss.wald.2 <- ceiling(temp$ss.eq2)

  # Compute sample sizes for score test
  temp <- getSampleSize(power, size.2sided, U20, var.score.null, var.score.alt)
  ss.score.1 <- ceiling(temp$ss.eq1)
  ss.score.2 <- ceiling(temp$ss.eq2)

  list(ss.wald.1=ss.wald.1, ss.wald.2=ss.wald.2, 
       ss.score.1=ss.score.1, ss.score.2=ss.score.2)

} # END: ss_data

######################################################################
# Main function for covariates, distribution
######################################################################
ss_covars_dist_f <- function(prev, logOR, rpdf=NULL, size.2sided=0.05, 
                  power=0.9, lam=0.5, interval=c(-7, 7), tol=0.0001,
                  rpdf.n=10000) {

  # rpdf   A Function to generate a random sample for the covariates 
  #          and exposure. Must be a vector of length beta

  if (is.null(rpdf)) {
    # Default to MVN(0, 1)
    rpdf <- function(n) {rmvnorm(n, mean=rep(0, length(logOR)))}
  }

  #  Check the pdf
  pdf.str <- check_pdf(rpdf)
  if (is.null(pdf.str)) pdf.str <- getObjName(rpdf)
 
  # Create new function called myrpdf and test it
  myrpdf <- NULL
  temp <- parsePDF(pdf.str, newPDF="myrpdf", newRandF="myNONFUNC") 
  eval(parse(text=temp$pdf))

  # Check rpdf
  vec <- myrpdf(1)
  if (length(vec) != length(logOR)) {
    stop("ERROR: length(logOR) must equal length of return vector from rpdf")
  }
 
  # Generate data
  mat <- myrpdf(rpdf.n)

  ret <- ss_data(prev, logOR, mat, size.2sided=size.2sided, power=power, 
                 lam=lam, interval=interval, tol=tol)

  ret
  
} # END: ss_covars_dist_f




