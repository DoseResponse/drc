\name{baro5}

\alias{baro5}

\title{The modified baro5 function}

\description{
  'baro5' allows specification of the baroreflex 5-parameter dose response function,
  under various constraints on the parameters.
}

\usage{
  baro5(fixed = c(NA, NA, NA, NA, NA), names = c("b1", "b2", "c", "d", "e"),
  method = c("1", "2", "3", "4"), ssfct = NULL)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b1, b2, c, d, e (see under 'Details').}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}                 
}

\details{
  The five-parameter function given by the expression

  \deqn{ y = c + \frac{d-c}{1+f\exp(b1(\log(x)-\log(e))) + (1-f)\exp(b2(\log(x)-\log(e)))}}

  \deqn{ f = 1/( 1 + \exp((2b1b2/|b1+b2|)(\log(x)-\log(e))))}
  
  If the difference between the parameters b1 and b2 is different from 0 then the function is asymmetric.
}

\value{
  The value returned is a list containing the nonlinear model function, the self starter function
  and the parameter names.    
}

\references{
  Ricketts, J. H. and Head, G. A. (1999) 
  A five-parameter logistic equation for investigating asymmetry of curvature in baroreflex studies. 
  \emph{Am. J. Physiol. (Regulatory Integrative Comp. Physiol. 46)}, \bold{277}, 441--454.
}

\author{Christian Ritz}

\note{See the example for the dataset \code{\link{heartrate}}.}

%\examples{}

\keyword{models}
\keyword{nonlinear}
