
###################################################################
# Main function for univariate case, ordinal exposure, distribution
###################################################################
pow_uni_dist_pmf <- function(prev, logOR, p, xvals, size.2sided=0.05, 
                     sampleSize=1000, lam=0.5, interval=c(-7, 7), tol=0.0001) {

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
  temp <- getPower(sampleSize, size.2sided, logOR, var.wald.null, var.wald.alt)
  pow.wald.1 <- temp$pow.eq1
  pow.wald.2 <- temp$pow.eq2

  # Compute sample sizes for score test
  mu1  <- getMean_pmf(xvals, px1) 
  mu0  <- getMean_pmf(xvals, px0) 
  temp <- getPower(sampleSize, size.2sided, mu1-mu0, var.score.null, var.score.alt)
  pow.score.1 <- temp$pow.eq1
  pow.score.2 <- temp$pow.eq2

  list(pow.wald.1=pow.wald.1, pow.wald.2=pow.wald.2, 
       pow.score.1=pow.score.1, pow.score.2=pow.score.2)

} # END: pow_uni_dist_pmf

