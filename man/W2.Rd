\name{W1.2}

\alias{W1.2}
\alias{w2}
\alias{W2.2}

\title{The two-parameter Weibull functions}

\description{
  'W1.2' is the two-parameter Weibull function where the lower limit is fixed at 0 and the upper limit
  is fixed at 1, mostly suitable for binomial/quantal responses.  
}

\usage{
  W1.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
  
  W2.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
}

\arguments{
  \item{upper}{numeric value. The fixed, upper limit in the model. Default is 1.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{additional arguments to be passed from the convenience functions.}  
}

\details{
  The two-parameter Weibull type 1 model is given by the expression
  \deqn{ f(x) = \exp(-\exp(b(\log(x)-\log(e)))).}
  
  The function is asymmetric about the inflection point, that is the parameter \eqn{\exp(e)}.
}

\value{
  See \code{\link{weibull1}}.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.
}

\seealso{
  Related functions are \code{\link{W1.3}}, \code{\link{W1.4}}, \code{\link{weibull1}} and \code{\link{weibull2}}.
}

\examples{

## Fitting a two-parameter Weibull model
earthworms.m1 <- drm(number/total~dose, weights = total,
data = earthworms, fct = W1.2(), type = "binomial")

summary(earthworms.m1)

}
\keyword{models}
\keyword{nonlinear}
