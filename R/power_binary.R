###################################################################
# Main function for univariate case, binary exposure, distribution
###################################################################
pow_uni_dist_bin <- function(prev, logOR, p, size.2sided=0.05, sampleSize=1000, lam=0.5,
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

  # Compute powers for Wald test
  temp <- getPower(sampleSize, size.2sided, logOR, var.wald.null, var.wald.alt)
  pow.wald.1 <- temp$pow.eq1
  pow.wald.2 <- temp$pow.eq2

  # Compute powers for score test
  temp <- getPower(sampleSize, size.2sided, px1-px0, var.score.null, var.score.alt)
  pow.score.1 <- temp$pow.eq1
  pow.score.2 <- temp$pow.eq2

  list(pow.wald.1=pow.wald.1, pow.wald.2=pow.wald.2, 
       pow.score.1=pow.score.1, pow.score.2=pow.score.2)

} # END: pow_uni_dist_bin

