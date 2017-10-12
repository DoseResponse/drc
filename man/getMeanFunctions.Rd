\name{getMeanFunctions}

\alias{getMeanFunctions}

\title{Display available dose-response models}

\description{
  Display information about available, built-in dose-response models.
}

\usage{
  getMeanFunctions(noParm = NA, fname = NULL, flist = NULL, display =TRUE)
}

\arguments{
  \item{noParm}{numeric specifying the number of parameters of the models to be displayed.
  The default (NA) results in display of all models, regardless of number of parameters.}
  \item{fname}{character string or vector of character strings specifying the short name(s) 
  of the models to be displayed (need to match exactly).}
  \item{flist}{list of built-in functions to be displayed.}
  \item{display}{logical indicating whether or not the requested models should be displayed on the R console.}  
}

\details{
  The arguments \code{noParm} and \code{fname} can be combined.
}

\value{
  An invisible list of functions or a list of strings with brief function descriptions is returned. 
}

%\references{}

\author{Christian Ritz}

%\note{}

%\seealso{}

\examples{

## Listing all functions
getMeanFunctions()

## Listing all functions with 4 parameters
getMeanFunctions(4)

## Listing all (log-)logistic functions
getMeanFunctions(fname = "L")

## Listing all three-parameter (log-)logistic or Weibull functions
getMeanFunctions(3, fname = c("LL", "W"))

## Listing all four-parameter (log-)logistic or Weibull functions
getMeanFunctions(4, fname = c("LL", "W"))

}
\keyword{models}
\keyword{nonlinear}
