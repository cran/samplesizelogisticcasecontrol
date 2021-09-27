% --- Source file: power_data.Rd ---
\name{power_data}
\alias{power_data}
\title{Power using pilot data}
\description{
Calculates the power of a case-control study with pilot data
}
\usage{power_data(prev, logOR, data, size.2sided=0.05, sampleSize=1000, cc.ratio=0.5,
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
               that gives a sample from the distribution of the confounders and exposure (see details).
              If a matrix or data frame, then the last column consists of random values for the exposure, 
              while the other columns are for the confounders. The order of the columns must match the order of 
              the vector \code{logOR}.
               The default is NULL.}
  \item{size.2sided}{Number between 0 and 1 giving the size of the 2-sided hypothesis test. The default is 0.05.}
  \item{sampleSize}{Sample size of the study (see details). The default is 1000.}
  \item{cc.ratio}{Number between 0 and 1 for the proportion of cases in the case-control sample.
                The default is 0.5.}
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

\details{
  The option \code{sampleSize} is not necessarily \code{nrow(data)}. The input \code{data} can be a
  small sample of pilot data that would be representative of the actual study data. 
}

\seealso{ \code{\link{power_binary}}, \code{\link{power_ordinal}}, \code{\link{power_continuous}} }

\examples{
  prev  <- 0.01
  logOR <- 0.3
  data  <- matrix(rnorm(100, mean=1.5), ncol=1)

  # Assuming exposuure is N(1.5, 1)
  power_data(prev, logOR, data) 

}
 