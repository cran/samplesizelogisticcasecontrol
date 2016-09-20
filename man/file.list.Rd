% --- Source file: file.list.Rd ---
\name{file.list}
\alias{file.list}
\docType{data}
\title{
 List to describe the covariate and exposure data 
}
\description{
 The list to describe the covariate and exposure data for the \code{data} option.
}

\format{
        The format is:
         List of 7
         \describe{
          \item{file}{ Data file containing the confounders and exposure variables.
                         No default.}
          \item{exposure}{ Name or column number in \code{file} for the exposure variable. 
                          This can also be a vector giving the columns to form an interaction
                          variable (see details).
                        No default. }
          \item{covars}{Character vector of variables names or numeric vector of column numbers
                        in \code{file} that will be confounders. These variables must be numeric.
                        The length and order of the \code{logOR} argument must match the length and 
                        order of c(\code{covars}, \code{exposure}).
                        The default is NULL.}
          \item{header}{ 0 or 1 if \code{file} has the first row as variable names. 
                         The default is determined from the first line of the \code{file}.}
          \item{delimiter}{ The delimiter in \code{file}.
                         The default is determined from the first two lines of the \code{file}.}
          \item{in.miss}{ Vector of character strings to define the missing values. This option
                          corresponds to the option \code{na.strings} in \code{read.table()}.
                          The default is "NA".}
          \item{subsetData}{ List of sublists to subset the data. 
                             Each sublist should contain the names "var", "operator" and "value" corresponding
                             to a variable name, operator and values of the variable.
                             Multiple sublists are logically connected by the AND operator. For example, \cr
                          subsetData=list(list(var="GENDER", operator="==", value="MALE")) \cr
                          will only include subjects with the string "MALE" for the GENDER variable. \cr
                          subsetData=list(list(var="AGE", operator=">", value=50), \cr
                                          list(var="STUDY", operator="\%in\%", value=c("A", "B", "C"))) \cr
                          will include subjects with AGE > 50 AND in STUDY A, B or C.  
                          The default is NULL.}
         }
}

\details{ In this list, \code{file} and \code{exposure} must be specified. If \code{exposure} is a vector
  of column names or column numbers, then an exposure variable will be created by multipling the columns
  defined in the vector to form the interaction variable. Thus, the columns must be numeric variables. 
  In this case, the length and order of \code{logOR} must match the length and order of
  c(\code{covars}, <new interaction variable>).
}
\keyword{misc}


