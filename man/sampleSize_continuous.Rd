% --- Source file: sampleSize_continuous.Rd ---
\name{sampleSize_continuous}
\alias{sampleSize_continuous}
\title{Sample size for a continuous exposure}
\description{
Calculates the required sample size of as case-control study with a continuous exposure variable
}
\usage{sampleSize_continuous(prev, logOR, distF=NULL, distF.support=c(-Inf, Inf), 
      data=NULL, size.2sided=0.05, power=0.9, cc.ratio=0.5, interval=c(-100, 100), 
      tol=0.0001, n.samples=10000, distF.var=NULL) 
}
\arguments{
  \item{prev}{Number between 0 and 1 giving the prevalence of disease. No default.}
  \item{logOR}{Vector of ordered log-odds ratios for the confounders and exposure.
               The last log-odds ratio in the vector is for the exposure. If the 
               option \code{data} (below) is specified, then the order must match the 
               order of \code{data}.  
               No default.}
  \item{distF}{NULL, a function or a character string giving the pdf of the exposure variable for the case 
               of no confounders, or giving the function to generate random vectors from the
               distribution formed by the confounders and exposure.
               For the case of no confounders, examples are \code{\link[stats]{dnorm}}, 
                "dnorm(x, mean=0.5, sd=2.1)", "dbeta(?, shape1=0.3, shape2=3)", "dchisq(whatever, df=1)".
                Notice that when \code{distF} is a character string, the first argument can be anything but
                must be given to serve as a place holder. 
               For the case of two confounders, an example might be a random vector generator from a 
                multivariate normal distribution "rmvnorm(x, c(0,0,0))". 
                User defined functions are also allowed, provided the user-defined function has only
                one input argument. The input argument would be a vector of quantiles if the user-defined 
                function is a pdf, or the input argument would be an integer specifiying the number
                of random vectors to generate if the user-defined function is a function to generate
                random vectors from the distribution of the confounders and exposure.
                An example pdf is the function H, where 
                  H <- function(x) \{ dunif(x, min=2, max=7) \}. 
                The default depends on other options (see details).}
  \item{distF.support}{Two element vector giving the domain of \code{distF}. This option is only used
                       when \code{distF} is a pdf.
                  The default is c(-Inf, Inf).}
  \item{data}{NULL, matrix, data frame or a list of type \code{\link{file.list}} 
               that gives a sample from the distribution of the confounders and exposure.
              If a matrix or data frame, then the last column consists of random values for the exposure, 
              while the other columns are for the confounders. The order of the columns must match the order of 
              the vector \code{logOR}.
               The default is NULL.}
  \item{size.2sided}{Number between 0 and 1 giving the size of the 2-sided hypothesis test. The default is 0.05.}
  \item{power}{Number between 0 and 1 for the desired power of the test. The default is 0.9.}
  \item{cc.ratio}{Number between 0 and 1 for the proportion of cases in the case-control sample. The default is 0.5.}
  \item{interval}{Two element vector giving the interval to search for the estimated intercept parameter. 
                  The default is c(-100, 100).}
  \item{tol}{Positive value giving the stopping tolerance for the root finding method to estimate
             the intercept parameter.
                  The default is 0.0001.}
  \item{n.samples}{Integer giving the number of random vectors to generate when the option \code{distF}
                  is specified and is a random vector generation function.
                  The default is 10000.}
  \item{distF.var}{The variance of the exposure variable for the case of no confounders. This option is 
                   for efficiency purposes. If not specified, the variance will be estimated by either
                   the empirical variance of a random sample from the distribution of the exposure or
                   by numerical integration.
                  The default is NULL.}
}
\value{
   A list containing four sample sizes, where two of them are for a Wald test and two for a score test.
   The two sample sizes for each test correspond to the equations for 
   \eqn{n_{1}}{n_{1}} and \eqn{n_{2}}{n_{2}}.
 }
\details{The \code{data} option takes precedance over the other options. If \code{data} is not specified,
  then the distribution of the exposure will be N(0,1) or MVN(0, 1) depending on whether there
  are confounders. 
}

\seealso{ \code{\link{sampleSize_binary}}, \code{\link{sampleSize_ordinal}}, \code{\link{sampleSize_data}} }

\examples{
  prev  <- 0.01
  logOR <- 0.3

  # No confounders, exposure assumed to be N(0,1)
  sampleSize_continuous(prev, logOR) 

  # Two confounders and exposure assumed to be MVN(0,1)
  beta <- c(0.1, 0.1, logOR)
  sampleSize_continuous(prev, beta) 

  # No confounders, exposure is beta(0.3, 3)
  sampleSize_continuous(prev, logOR, distF="dbeta(m, shape1=0.3, shape2=3)",
                        distF.support=c(0, 1)) 
}
 
