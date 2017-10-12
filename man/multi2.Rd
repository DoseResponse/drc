\name{multi2}

\alias{multi2}

\title{
  Multistage dose-response model with quadratic terms
}

\description{
   The multistage dose-response model is a combination of log-logistic models that should be useful for describing
   more complex dose-response patterns.
}

\usage{
  multi2(
  fixed = c(NA, NA, NA, NA, NA), 
  names = c("b1", "b2", "b3", "c", "d"), 
  ssfct = NULL, 
  fctName, 
  fctText)
}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at what value they are fixed. 
  NAs are used for parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage').}
  \item{ssfct}{a self starter function to be used.}    
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}
}

\details{
  The multistage model function with quadratic terms is defined as follows
  
  \deqn{ f(x) = c + (d-c)\exp(-b1-b2x-b3x^2)}
  
  where x denotes the dose or the logarithm-transformed dose.
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}

\references{
  Wheeler, M. W., Bailer, A. J. (2009)
  Comparing model averaging with other model selection strategies for benchmark dose estimation,
  \emph{Environmental and Ecological Statistics}, \bold{16}, 37--51.
}

\author{
  Christian Ritz
}

%\note{}

%\seealso{}

%\examples{}

\keyword{models}
\keyword{nonlinear}
