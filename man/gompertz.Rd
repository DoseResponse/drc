\name{gompertz}

\alias{gompertz}
\alias{G.2}
\alias{G.3}
\alias{G.3u}
\alias{G.4}

\title{Mean function for the Gompertz dose-response or growth curve}

\description{
  This function provides a very general way of specifying the mean function of the decreasing or incresing 
  Gompertz dose-response or growth curve models.
}

\usage{
  gompertz(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
  method = c("1", "2", "3", "4"), ssfct = NULL,
  fctName, fctText)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The order of the parameters is: b, c, d, e (see under 'Details' for the precise meaning of each parameter).}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{character string used internally by convenience functions (optional).}  
  \item{fctText}{character string used internally by convenience functions (optional).}  
}

\details{
  The Gompertz model is given by the mean function
  
  \deqn{ f(x) = c + (d-c)(\exp(-\exp(b(x-e)))) }  

  and it is a dose-response/growth curve on the entire real axis, that is it is not limited to
  non-negative values even though this is the range for most dose-response and growth data. One consequence is
  that the curve needs not reach the lower asymptote at dose 0.
  
  If \deqn{b<0} the mean function is increasing and it is decreasing for \deqn{b>0}. The decreasing Gompertz model is
  not a well-defined dose-response model and other dose-response models such as the Weibull models 
  should be used instead.
    
  Various re-parameterisations of the model are used in practice.
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

\references{
  Seber, G. A. F. and Wild, C. J. (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (p. 331).
}

\author{Christian Ritz}

\note{
  The function is for use with the function \code{\link{drm}}, but typically the convenience functions
  \code{\link{G.2}}, \code{\link{G.3}}, \code{\link{G.3u}}, and \code{\link{G.4}} should be used.
}

\seealso{The Weibull model \code{\link{weibull2}} is closely related to the Gompertz model.}

%\examples{}

\keyword{models}
\keyword{nonlinear}

\concept{Gompertz}