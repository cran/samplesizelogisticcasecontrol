% --- Source file: sampleSize.Rd ---
\name{samplesizelogisticcasecontrol-package}
\alias{samplesizelogisticcasecontrol}
\docType{package}
\title{
Sample size and power calculations for case-control studies
}
\description{This package can be used to calculate the required sample size needed for 
  case-control studies to have sufficient power, or calculate the power of a
  case-control study for a given sample size. 
  To calculate the sample size, one needs to
  specify the significance level \eqn{\alpha}{alpha}, power \eqn{\gamma}{gamma}, and the 
  hypothesized non-null \eqn{\theta}{theta}. Here \eqn{\theta}{theta} is a log odds ratio
  for an exposure main effect or \eqn{\theta}{theta} is an interaction effect on the 
  logistic scale.
  Choosing \eqn{\theta}{theta} requires subject matter knowledge to understand how
  strong the association needs to be to have practical importance.   
Sample size varies inversely with
  \eqn{\theta^{2}}{theta^{2}} and is thus highly dependent on \eqn{\theta}{theta}. 
}
\details{
The main functions in the package are for different types of exposure variables, where the
exposure variable is the variable of interest in a hypothesis test. 
The functions \code{\link{sampleSize_binary}} and \code{\link{power_binary}}
 can be used for a binary exposure variable (X = 0 or 1), 
while the functions \code{\link{sampleSize_ordinal}} and \code{\link{power_ordinal}} 
 is a more general function that can be used for 
an ordinal exposure variable (X takes the values 0, 1, ..., k). 
\code{\link{sampleSize_continuous}} and \code{\link{power_continuous}} are useful for a continuous exposure variable and
\code{\link{sampleSize_data}} and \code{\link{power_data}} can be used when pilot data is available that defines
the distribution of the exposure and other confounding variables. Each function will return the
sample sizes or power for a Wald-type test and a score test. When there are no adjustments for confounders, the user can 
specify a general distribution for the exposure variable. With confounders, either pilot data or a function to 
generate random samples from the multivariate distribution of the confounders and exposure variable must
be given. 


If the parameter of interest, \eqn{\theta}{theta},
  is one dimensional, then the test statistic is often asymptotically equivalent
  to a test of the form 
   \eqn{T > Z_{1-\alpha}\sigma_{0}n^{-\frac{1}{2}}}{T > Z_{1-alpha}sigma_{0}n^{-1/2}} or
   \eqn{T > Z_{1-\alpha}\sigma_{\theta}n^{-\frac{1}{2}}}{T > Z_{1-alpha}sigma_{theta}n^{-1/2}}, where
  \eqn{Z_{1-\alpha}}{Z_{1-alpha}} is the \eqn{1-\alpha}{1-alpha} quantile of a standard
  normal distribution, \eqn{n}{n} is the total sample size (cases plus controls), 
   and \eqn{n^{\frac{1}{2}}T}{n^{1/2}T} is 
  normally distributed with mean 0 and null variance \eqn{\sigma_{0}^{2}}{sigma_{0}^{2}}.
  Depending on which critical value
  \eqn{Z_{1-\alpha}\sigma_{0}n^{-\frac{1}{2}}}{Z_{1-alpha}sigma_{0}n^{-1/2}} or
  \eqn{Z_{1-\alpha}\sigma_{\theta}n^{-\frac{1}{2}}}{Z_{1-alpha}sigma_{theta}n^{-1/2}} 
 of the test was used, the formulas for sample size are obtained by inverting the
  equations for power:

  \eqn{n_{1} = (Z_{\gamma}\sigma_{\theta} + Z_{1-\alpha}\sigma_{0})^{2}/\theta^{2}}{n_{1} = (Z_{gamma}sigma_{theta} + Z_{1-alpha}sigma_{0})^{2}/theta^{2}}
 or
\eqn{n_{2} = (Z_{\gamma} + Z_{1-\alpha})^{2}\sigma_{\theta}^{2}/\theta^{2}}{n_{2} = (Z_{gamma} + Z_{1-alpha})^{2}sigma_{theta}^{2}/theta^{2}}. 


}
\author{Mitchell H. Gail <gailm@mail.nih.gov>}
\keyword{package}
\references{
Gail, M.H. and Haneuse, S. Power and sample size for case-control studies.  
In Handbook of Statistical Methods for Case-Control Studies.  
Editors:  Ornulf Borgan, Norman Breslow, Nilanjan Chatterjee, Mitchell Gail, Alastair Scott, Christopher Wild.  
Chapman and Hall/CRC, Taylor and Francis Group, New York, 2018, pages 163-187.

Gail, M. H and Haneuse, S. Power and sample size for multivariate logistic modeling of unmatched case-control studies. 
Statistical Methods in Medical Research. 2019;28(3):822-834, \cr
 https://doi.org/10.1177/0962280217737157
}

