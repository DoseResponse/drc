\name{cedergreen}

\alias{cedergreen}
\alias{CRS.6}
\alias{ucedergreen}

\title{The Cedergreen-Ritz-Streibig model}

\description{
  'cedergreen' provides a very general way of specifying then Cedergreen-Ritz-Streibig 
  modified log-logistic model for describing hormesis, under various constraints on the parameters.

  \code{\link{CRS.6}} is the extension of \code{link{cedergreen}} with freely varying alpha parameter.

  For u-shaped hormesis data 'ucedergreen' provides a very general way of specifying the 
  Cedergreen-Ritz-Streibig modified log-logistic model, under various constraints on the parameters.
}

\usage{
  cedergreen(fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), 
  method = c("1", "2", "3", "4"), ssfct = NULL,
  alpha, fctName, fctText)
  
  CRS.6(fixed = c(NA, NA, NA, NA, NA, NA), 
  names = c("b","c","d","e","f","g"),
  method = c("1", "2", "3", "4"), ssfct = NULL)
  
  ucedergreen(fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"), 
  method = c("1", "2", "3", "4"), ssfct = NULL,
  alpha)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}                 
  \item{alpha}{numeric value between 0 and 1, reflecting the steepness of the hormesis peak. 
  This argument needs to be specified.}
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}  
}

\details{
The model is given by the expression
  
  \deqn{ f(x) = c + \frac{d-c+f exp(-1/(x^{\alpha}))}{1+exp(b(log(x)-log(e)))}}
  
which is a five-parameter model (alpha is fixed or freely varying). Not all features (eg EC/ED calculation) 
are available for the model with freely varying alpha.
  
It is a modification of the four-parameter logistic curve to take hormesis into account.
  
The u-shaped model is given by the expression
  
  \deqn{ f(x) = cd - \frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

\references{
  Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) 
  Improved empirical models describing hormesis, 
  \emph{Environmental Toxicology and Chemistry} 
  \bold{24}, 3166--3172.
}

\author{Christian Ritz}

\note{
  The functions are for use with the functions \code{\link{drm}}.
}

\seealso{
  For fixed alpha, several special cases are handled by the following convenience functions 
  \code{\link{CRS.4a}}, \code{\link{CRS.4b}}, 
  \code{\link{CRS.4c}}, \code{\link{CRS.5a}}, \code{\link{CRS.5b}}, \code{\link{CRS.5c}},
  \code{\link{UCRS.4a}}, \code{\link{UCRS.4b}}, \code{\link{UCRS.4c}}, \code{\link{UCRS.5a}},
  \code{\link{UCRS.5b}}, \code{\link{UCRS.5c}} where a, b and c correspond to 
  the pre-specified alpha values 1, 0.5 and 0.25, respectively.
}

\examples{

## Estimating CRS model with alpha unknown
lettuce.crsm1 <- drm(weight~conc, data = lettuce, fct = CRS.6())
summary(lettuce.crsm1)
plot(lettuce.crsm1)  # oops: not increasing until hormesis peak

}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation u-shaped}
