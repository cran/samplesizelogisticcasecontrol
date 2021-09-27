
###################################################################
# Main function for pilot data, covariates
###################################################################
pow_data <- function(prev, beta, data, size.2sided=0.05, sampleSize=1000, lam=0.5,
                    interval=c(-7, 7), tol=0.0001) {

  nbeta      <- length(beta)
  ncov       <- nbeta - 1
  data       <- as.matrix(data)
 
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

  # Compute powers for Wald test
  temp <- getPower(sampleSize, size.2sided, beta[nbeta], var.wald.null, var.wald.alt)
  pow.wald.1 <- temp$pow.eq1
  pow.wald.2 <- temp$pow.eq2

  # Compute powers for score test
  temp <- getPower(sampleSize, size.2sided, U20, var.score.null, var.score.alt)
  pow.score.1 <- temp$pow.eq1
  pow.score.2 <- temp$pow.eq2

  list(pow.wald.1=pow.wald.1, pow.wald.2=pow.wald.2, 
       pow.score.1=pow.score.1, pow.score.2=pow.score.2)

} # END: pow_data

######################################################################
# Main function for covariates, distribution
######################################################################
pow_covars_dist_f <- function(prev, logOR, rpdf=NULL, size.2sided=0.05, 
                  sampleSize=1000, lam=0.5, interval=c(-7, 7), tol=0.0001) {

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
  mat <- myrpdf(sampleSize)

  ret <- pow_data(prev, logOR, mat, size.2sided=size.2sided, 
                 lam=lam, interval=interval, tol=tol)

  ret
  
} # END: pow_covars_dist_f




