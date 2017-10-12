\name{gammadr}

\Rdversion{1.1}

\alias{gammadr}

\title{
  Gamma dose-response model
}

\description{
   The gamma dose-response model is a four-parameter model derived from the cumulative distribution function of the gamma distribution.
}

\usage{
  gammadr(fixed = c(NA, NA, NA, NA), 
  names = c("b", "c", "d", "e"), fctName, fctText)
}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at what value they are fixed. 
  NAs are used for parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage').}
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}
}

\details{
  Following Wheeler and Bailer (2009) the model function is defined as follows:
  
  \deqn{f(x) = c + (d-c) * pgamma(b*x, e, 1)}
  
  
  This model is only suitable for increasing dose-response data.
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

%\seealso{
%  The basic component in the two-phase model is the log-logistic model
%  \code{\link{llogistic}}.
%}

%\examples{}

\keyword{models}
\keyword{nonlinear}
