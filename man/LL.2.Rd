\name{LL.2}

\alias{LL.2}
\alias{l2}

\alias{LL2.2}

\title{The two-parameter log-logistic function}

\description{
  'LL.2' and 'LL2.2' provide the two-parameter log-logistic function where the lower limit is fixed at 0 and the upper limit
  is fixed at 1, mostly suitable for binomial/quantal responses. 
}

\usage{
  LL.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
  
  l2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
  
  LL2.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)  
}

\arguments{
  \item{upper}{numeric value. The fixed, upper limit in the model. Default is 1.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters. The default is reasonable.}
  \item{...}{Additional arguments (see \code{\link{llogistic}}).}
}

\details{
  The two-parameter log-logistic function is given by the expression
  \deqn{ f(x) = \frac{1}{1+\exp(b(\log(x)-\log(e)))}}
  
  or in another parameterisation
  \deqn{ f(x) = \frac{1}{1+\exp(b(\log(x)-e))}}  
  
  The model function is symmetric about the inflection point (\eqn{e}).  
}

\value{
  See \code{\link{llogistic}}.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.
}

\seealso{
  Related functions are \code{\link{LL.3}}, \code{\link{LL.4}}, \code{\link{LL.5}} and the more general 
  \code{\link{llogistic}}.
}

\examples{

## Fitting a two-parameter logistic model 
##  to binomial responses (a logit model)
earthworms.m1 <- drm(number/total~dose, weights=total, 
data = earthworms, fct = LL.2(), type = "binomial")

plot(earthworms.m1)  # not fitting at the upper limit!

}

\keyword{models}
\keyword{nonlinear}

\concept{Hill logit 2-parameter two-parameter}