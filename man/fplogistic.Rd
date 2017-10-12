\name{fplogistic}

\Rdversion{1.1}

\alias{fplogistic}
\alias{FPL.4}

\title{Fractional polynomial-logistic dose-response models}

\description{
  Model function for specifying dose-response models that are a combination of a logistic model and an appropriate 
  class of fractional polynomials.
}

\usage{
fplogistic(p1, p2, fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
method = c("1", "2", "3", "4"), ssfct = NULL, fctName, fctText)

FPL.4(p1, p2, fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
}

\arguments{
  \item{p1}{numeric denoting the negative power of log(dose+1) in the fractional polynomial.}
  \item{p2}{numeric denoting the positive power of log(dose+1) in the fractional polynomial.}
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
    NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f 
    (see under 'Details').}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.} 
  \item{...}{Additional arguments (see \code{\link{fplogistic}}).} 
}

\details{
  The fractional polynomial dose-response models introduced by Namata \emph{et al.} (2008) are implemented
  using the logistic model as base.
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}

\references{
  Namata, Harriet and Aerts, Marc and Faes, Christel and Teunis, Peter (2008)
  Model Averaging in Microbial Risk Assessment Using Fractional Polynomials,
  \emph{Risk Analysis} \bold{28}, 891--905.
}

\author{Christian Ritz}

%\note{}

\seealso{
  Examples are found on the hep page of \code{\link{maED}}.
}

%\examples{}

\keyword{models}
\keyword{nonlinear}