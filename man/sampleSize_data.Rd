% --- Source file: sampleSize_data.Rd ---
\name{sampleSize_data}
\alias{sampleSize_data}
\title{Sample size using pilot data}
\description{
Calculates the required sample size of a case-control study with pilot data
}
\usage{sampleSize_data(prev, logOR, data, size.2sided=0.05, power=0.9, cc.ratio=0.5,
        interval=c(-100, 100), tol=0.0001) 
}
\arguments{
  \item{prev}{Number between 0 and 1 giving the prevalence of disease. No default.}
  \item{logOR}{Vector of ordered log-odds ratios for the confounders and exposure.
               The last log-odds ratio in the vector is for the exposure. If the 
               option \code{data} (below) is specified, then the order must match the 
               order of \code{data}.  
               No default.}
  \item{data}{Matrix, data frame or a list of type \code{\link{file.list}} 
               that gives a sample from the distribution of the confounders and exposure.
              If a matrix or data frame, then the last column consists of random values for the exposure, 
              while the other columns are for the confounders. The order of the columns must match the order of 
              the vector \code{logOR}.
               The default is NULL.}
  \item{size.2sided}{Number between 0 and 1 giving the size of the 2-sided hypothesis test. The default is 0.05.}
  \item{power}{Number between 0 and 1 for the desired power of the test. The default is 0.9.}
  \item{cc.ratio}{Number between 0 and 1 for the proportion of cases in the case-control sample.
                The default is 0.5.}
  \item{interval}{Two element vector giving the interval to search for the estimated intercept parameter. 
                  The default is c(-100, 100).}
  \item{tol}{Positive value giving the stopping tolerance for the root finding method to estimate
             the intercept parameter.
                  The default is 0.0001.}
}
\value{
   A list containing four sample sizes, where two of them are for a Wald test and two for a score test.
   The two sample sizes for each test correspond to the equations for 
   \eqn{n_{1}}{n_{1}} and \eqn{n_{2}}{n_{2}}.
 }

\seealso{ \code{\link{sampleSize_binary}}, \code{\link{sampleSize_ordinal}}, \code{\link{sampleSize_continuous}} }

\examples{
  prev  <- 0.01
  logOR <- 0.3
  data  <- matrix(rnorm(100, mean=1.5), ncol=1)

  # Assuming exposuure is N(1.5, 1)
  sampleSize_data(prev, logOR, data) 

}
 