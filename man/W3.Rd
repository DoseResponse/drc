\name{W1.3}

\alias{W1.3}
\alias{w3}
\alias{W2.3}
\alias{W2x.3}

\alias{W1.3u}
\alias{W2.3u}

\title{The three-parameter Weibull functions}

\description{
  'W1.3' and \code{W2.3} provide the three-parameter Weibull function, self starter function and names of the parameters.
  
  'W1.3u' and 'W2.3u' provide three-parameter Weibull function where the upper limit is equal to 1, mainly
  for use with binomial/quantal response.     
}

\usage{
  W1.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  W2.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  W2x.3(fixed = c(NA, NA, NA), names = c("d", "e", "t0"), ...)
  
  W1.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)  
  
  W2.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)    
}

\arguments{
  \item{upper}{numeric value. The fixed, upper limit in the model. Default is 1.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{additional arguments to be passed from the convenience functions.}   
}

\details{
  The three-parameter Weibull type 1 model is given by the expression
  \deqn{ f(x) = 0 + (d-0)\exp(-\exp(b(\log(x)-\log(e)))).}
  
  The model function is asymmetric about the inflection point, which is the parameter \eqn{\exp(e)}.

  The three-parameter Weibull type 1 model with upper limit 1 is given by the expression
  \deqn{ f(x) = 0 + (1-0)\exp(-\exp(b(\log(x)-\log(e)))).}

}

\value{
  See \code{\link{weibull1}}.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{This function is for use with the function \code{\link{drm}}.}

\seealso{Related functions are \code{\link{W1.4}} and \code{\link{weibull1}}.}

\examples{

## Fitting a three-parameter Weibull model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
ryegrass.m1

}
\keyword{models}
\keyword{nonlinear}
