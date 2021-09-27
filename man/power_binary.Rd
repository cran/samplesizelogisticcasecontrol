% --- Source file: power_binary.Rd ---
\name{power_binary}
\alias{power_binary}
\title{Power for a binary exposure}
\description{
Calculates the power of as case-control study with a binary exposure variable
}
\usage{
power_binary(prev, logOR, probXeq1=NULL, distF=NULL, data=NULL, 
      size.2sided=0.05, sampleSize=1000, cc.ratio=0.5, interval=c(-100, 100), tol=0.0001) 
}
\arguments{
  \item{prev}{Number between 0 and 1 giving the prevalence of disease. No default.}
  \item{logOR}{Vector of ordered log-odds ratios for the confounders and exposure.
               The last log-odds ratio in the vector is for the exposure. If the 
               option \code{data} (below) is specified, then the order must match the 
               order of \code{data}.  
               No default.}
  \item{probXeq1}{NULL or a number between 0 and 1 giving the probability that the exposure
                  variable is 1. If set to NULL, the the \code{data} option must be specified so 
                 that \code{probXeq1} can be estimated. The default is NULL.}
  \item{distF}{NULL, a function or a character string giving the function to generate random
                vectors from the distribution of the confounders and exposure. The order of the returned
                vector must match the order of \code{logOR}.
                User defined functions are also allowed, provided the user-defined function has only
                one integer valued argument that inputs the number of random vectors to generate. 
                For instance the header of a user-defined function called "userF" would be 
                userF <- function(n). 
                The default depends on other options (see details).}
  \item{data}{NULL, matrix, data frame or a list of type \code{\link{file.list}} that gives
       a sample from the distribution of the confounders and exposure.
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
\details{If there are no confounders (length(logOR) = 1), then either \code{probXeq1} or \code{data} must
  be specified, where \code{probXeq1} takes precedance. If there are confounders (length(logOR) > 1), then
  either \code{data} or \code{distF} must be specified, where \code{data} takes precedance.
}

\seealso{ \code{\link{power_continuous}}, \code{\link{power_ordinal}}, \code{\link{power_data}} }

\examples{
  prev  <- 0.01
  logOR <- 0.3

  # No confounders, Prob(X=1)=0.2
  power_binary(prev, logOR, probXeq1=0.2) 

  # Generate data for a N(0,1) confounder and binary exposure
  data <- cbind(rnorm(1000), rbinom(1000, 1, 0.4))
  beta <- c(0.1, 0.2)
  power_binary(prev, beta, data=data) 

  # Define a function to generate random vectors for two confounders and the binary exposure
  f <- function(n) {cbind(rnorm(n), rbinom(n, 3, 0.5), rbinom(n, 1, 0.3))}
  logOR <- c(0.2, 0.3, 0.25)
  power_binary(prev, logOR, distF=f) 

}
 

\keyword{ power } 
