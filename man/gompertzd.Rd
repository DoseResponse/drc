\name{gompertzd}

\alias{gompertzd}

\title{The derivative of the Gompertz function}

\description{
  'gompertzd' provides a way of specifying the derivative of the Gompertz function
  as a dose-response model.
}

\usage{
  gompertzd(fixed = c(NA, NA), names = c("a", "b"))
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is (notice the order): a, b.}
}

\details{
  The derivative of the Gompertz function is defined as
  
  \deqn{ f(x) = a \exp(bx-a/b(exp(bx)-1))}
  
  For \eqn{a>0} and \eqn{b} not 0, the function is decreasing, equaling \eqn{a} at \eqn{x=0} 
  and approaching 0 at plus infinity.
}

\value{
  The value returned is a list containing the model function, the self starter function
  and the parameter names.
}

%\references{}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.
}

%\seealso{}

%\examples{}

\keyword{models}
\keyword{nonlinear}
