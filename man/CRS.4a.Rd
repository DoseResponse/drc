\name{CRS.4a}

\alias{CRS.4a}
\alias{CRS.4b}
\alias{CRS.4c}

\alias{ml3a}
\alias{ml3b}
\alias{ml3c}

\alias{UCRS.4a}
\alias{UCRS.4b}
\alias{UCRS.4c}

\alias{uml3a}
\alias{uml3b}
\alias{uml3c}

\title{The Cedergreen-Ritz-Streibig model}

\description{
  'CRS.4a', 'CRS.4b' and 'CRS.4c' provide the Cedergreen-Ritz-Streibig modified log-logistic model 
  for describing hormesis with the lower limit equal to 0.

  'UCRS.4a', 'UCRS.4b' and 'UCRS.4c' provide the Cedergreen-Ritz-Streibig modified log-logistic model 
  for describing u-shaped hormesis with the lower limit equal to 0.
}

\usage{
  CRS.4a(names = c("b", "d", "e", "f"), ...)

  UCRS.4a(names = c("b", "d", "e", "f"), ...)
}

\arguments{
  \item{names}{a vector of character strings giving the names of the parameters. 
  The default is reasonable (see above).}
  \item{...}{additional arguments to be passed from the convenience functions.}
}

\details{
  The model is given by the expression
  \deqn{ f(x) = 0 + \frac{d-0+f \exp(-1/x)}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account.

  The u-shaped model is given by the expression
  \deqn{ f(x) = 0 + d - \frac{d-0+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
  The a,b,c models are obtained by setting alpha equal to 1, 0.5 and 0.25, respectively.    
}

\value{
  See \code{\link{cedergreen}}.
}

\references{
  See the reference under \code{\link{cedergreen}}.
}

\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.
}

\seealso{
  Similar functions are \code{\link{CRS.5a}} and \code{\link{UCRS.5a}}, 
  but with an extra parameter for the lower limit.}

\examples{

## Fitting modified logistic models
lettuce.crsm1 <- drm(lettuce[,c(2,1)], fct=CRS.4a())
summary(lettuce.crsm1)
ED(lettuce.crsm1, c(50))

## Need to explicitly specify that the upper limit
## is the reference in order to get ED10 and ED90 right
ED(lettuce.crsm1, c(10, 50, 90), reference = "upper")

lettuce.crsm2 <- drm(lettuce[,c(2,1)], fct=CRS.4b())
summary(lettuce.crsm2)
ED(lettuce.crsm2, c(50))

lettuce.crsm3 <- drm(lettuce[,c(2,1)], fct=CRS.4c())
summary(lettuce.crsm3)
ED(lettuce.crsm3, c(50))

}

\keyword{models}
\keyword{nonlinear}
