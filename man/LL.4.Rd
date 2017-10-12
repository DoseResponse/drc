\name{LL.4}

\alias{LL.4}
\alias{l4}

\alias{LL2.4}

\title{The four-parameter log-logistic function}

\description{
  'LL.4' and 'LL2.4' provide the four-parameter log-logistic function, self starter function, names of the parameters and, optionally,
  first and second derivatives for a faster estimation.  
}

\usage{
  LL.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)

  l4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
  
  LL2.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)  
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{Additional arguments (see \code{\link{llogistic}}).}
}

\details{
  The four-parameter log-logistic function is given by the expression
  \deqn{ f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}}
  
  or in another parameterisation (converting the term \eqn{\log(e)} into a parameter)
  \deqn{ f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\tilde{e}))}}  
  
  The function is symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{llogistic}}.
}

\references{
   Seber, G. A. F. and Wild, C. J (1989) \emph{Nonlinear Regression}, New York: Wiley \& Sons (p. 330).
}

\author{Christian Ritz and Jens C. Streibig}

\note{This function is for use with the function \code{\link{drm}}.}

\seealso{Setting \eqn{c=0} yields \code{\link{LL.3}}. See also \code{\link{LL.5}}.}

\examples{

spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
spinach.m1


}

\keyword{models}
\keyword{nonlinear}

\concept{Hill 4-parameter four-parameter}
