% --- Source file: power_ordinal.Rd ---
\name{power_ordinal}
\alias{power_ordinal}
\title{Power for an ordinal exposure}
\description{
Calculates the power of as case-control study with an ordinal exposure variable
}
\usage{power_ordinal(prev, logOR, probX=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, sampleSize=1000, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) 
}
\arguments{
  \item{prev}{Number between 0 and 1 giving the prevalence of disease. No default.}
  \item{logOR}{Vector of ordered log-odds ratios per category increase for the confounders and exposure.
               The last log-odds ratio in the vector is for the exposure. 
               If the exposures are coded 0, 1, ..., k (k+1 categories), then the \code{logOR}
               corresponds to a one category increase.
               If the option \code{data} (below) is specified, then the order must match the 
               order of \code{data}.  
               No default.}
  \item{probX}{NULL or a vector that sums to 1 giving the probability that the exposure
                  variable is equal to i, i = 0, 1, ..., k. 
                  If set to NULL, the the \code{data} option must be specified so 
                 that \code{probX} can be estimated. The default is NULL.}
  \item{distF}{NULL, a function or a character string giving the function to generate random
                vectors from the distribution of the confounders and exposure. The order of the returned
                vector must match the order of \code{logOR}.
                The default depends on other options (see details).}
  \item{data}{NULL, matrix, data frame or a list of type \code{\link{file.list}} 
             that gives a sample from the distribution of the confounders and exposure.
              If a matrix or data frame, then the last column consists of random values for the exposure, 
              while the other columns are for the confounders. The order of the columns must match the order of 
              the vector \code{logOR}.
               The default is NULL.}
  \item{size.2sided}{Number between 0 and 1 giving the size of the 2-sided hypothesis test. The default is 0.05.}
  \item{sampleSize}{Sample size of the study. The default is 1000.}
  \item{cc.ratio}{Number between 0 and 1 for the proportion of cases in the case-control sample. The default is 0.5.}
  \item{interval}{Two element vector giving the interval to search for the estimated intercept parameter. 
                  The default is c(-100, 100).}
  \item{tol}{Positive value giving the stopping tolerance for the root finding method to estimate
             the intercept parameter.
                  The default is 0.0001.}
}
\value{
   A list containing four powers, where two of them are for a Wald test and two for a score test.
   The two powers for each test correspond to the equations for 
   \eqn{n_{1}}{n_{1}} and \eqn{n_{2}}{n_{2}}.
 }
\details{If there are no confounders (length(logOR) = 1), then either \code{probX} or \code{data} must
  be specified, where \code{probX} takes precedance. If there are confounders (length(logOR) > 1), then
  either \code{data} or \code{distF} must be specified, where \code{data} takes precedance.
}

\seealso{ \code{\link{power_continuous}}, \code{\link{power_binary}}, \code{\link{power_data}} }

\examples{
  prev  <- 0.01
  logOR <- 0.3

  # No confounders, Prob(X=1)=0.2
  power_ordinal(prev, logOR, probX=c(0.8, 0.2)) 

  # Generate data for a N(0,1) confounder and ordinal exposure with 3 levels
  data <- cbind(rnorm(1000), rbinom(1000, 2, 0.5))
  beta <- c(0.1, 0.2)
  power_ordinal(prev, beta, data=data) 

  # Define a function to generate random vectors for two confounders and an ordinal
  #   exposure with 5 levels
  f <- function(n) {cbind(rnorm(n), rbinom(n, 1, 0.5), rbinom(n, 4, 0.5))}
  beta <- c(0.2, 0.3, 0.25)
  power_ordinal(prev, beta, distF=f) 

}
 
