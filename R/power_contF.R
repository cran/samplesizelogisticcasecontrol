
######################################################################
# Main function for univariate case, continuous exposure, distribution
######################################################################
pow_uni_dist_f <- function(prev, logOR, pdf, pdf.range, size.2sided=0.05, 
                  sampleSize=1000, lam=0.5, interval=c(-7, 7), tol=0.0001,
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
  #varX     <- getVarX(pdf.var, rpdf.str, mypdf, n=pdf.var.n, lower=pdf.range[1],
  #                    upper=pdf.range[2]) 
  beta     <- logOR

  # Get the intercept in general population
  alphastar <- getPopInt_dist_f(beta, prev, mypdf, pdf.range, interval=interval, tol=tol) 

  # Function to compute the intercept in case-control population
  alphacc <- getCCInt_dist_f(alphastar, beta, prev, lam, mypdf, pdf.range, f3Cont, 
                     interval=interval, tol=tol)

  # Compute information matrix
  info.alt <- getInfo_dist_f_alt(alphacc, alphastar, beta, prev, lam, pdf.range, mypdf)

  a    <- pdf.range[1]
  b    <- pdf.range[2]
  mux1 <- myintegrate_PDF(f5,a,b,method=2,alphastar,beta,mypdf)/prev # mean x in cases
  mux0 <- myintegrate_PDF(f6,a,b,method=2,alphastar,beta,mypdf)/(1-prev) #mean x in controls

  # For the NULL
  mu1out <- myintegrate_PDF(mux, a,b,method=2,alphacc,alphastar,beta,prev,lam,mypdf)
  mu2out <- myintegrate_PDF(mux2,a,b,method=2,alphacc,alphastar,beta,prev,lam,mypdf)

  # Get the variances for the Wald and score test
  temp           <- getVar_dist_f(lam, info.alt, mu1out, mu2out)
  var.wald.null  <- temp$var.wald.null
  var.wald.alt   <- temp$var.wald.alt
  var.score.null <- temp$var.score.null
  var.score.alt  <- temp$var.score.alt

  # Compute sample sizes for Wald test
  temp <- getPower(sampleSize, size.2sided, logOR, var.wald.null, var.wald.alt)
  pow.wald.1 <- temp$pow.eq1
  pow.wald.2 <- temp$pow.eq2

  # Compute sample sizes for score test
  temp <- getPower(sampleSize, size.2sided, mux1-mux0, var.score.null, var.score.alt)
  pow.score.1 <- temp$pow.eq1
  pow.score.2 <- temp$pow.eq2

  list(pow.wald.1=pow.wald.1, pow.wald.2=pow.wald.2, 
       pow.score.1=pow.score.1, pow.score.2=pow.score.2)

} # END: pow_uni_dist_f

