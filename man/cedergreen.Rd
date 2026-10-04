\name{cedergreen}

\alias{cedergreen}
\alias{CRS.6}
\alias{ucedergreen}

\title{The Cedergreen-Ritz-Streibig model}

\description{
  'cedergreen' provides a very general way of specifying the Cedergreen-Ritz-Streibig 
  modified log-logistic model for describing inverse U-shaped and U-shaped hormesis, 
  under various constraints on the parameters.

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
  alpha, fctName, fctText)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}                 
  \item{alpha}{numeric value between 0 and 1, reflecting the steepness of the hormesis peak.}
  \item{fctName}{optional character string, used internally for naming convenience functions.}  
  \item{fctText}{optional character string, used internally for naming convenience functions.}  
}

\details{
The inverse U-shaped model is given by the expression
  
  \deqn{ f(x) = c + \frac{d-c+f exp(-1/(x^{\alpha}))}{1+exp(b(log(x)-log(e)))}}
  
which is a five-parameter model (alpha is fixed or freely varying). Not all features (eg EC/ED calculation) 
are available for the model with freely varying alpha.
  
It is a modification of the four-parameter log-logistic model function to take hormesis into account.
  
The corresponding U-shaped model is given by the expression
  
  \deqn{ f(x) = d - \frac{d-c + f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}}
  
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
  https://doi.org/10.1897/05-014R.1

  Sweeney, C. J. and Bottoms, M. and Hutcheson, K. and Schulz, L. (2026) 
  Dose-response OECD 216 studies for agrochemicals, 
  \emph{Integrated Environmental Assessment and Management}, vjag057. 
  https://doi.org/10.1093/inteam/vjag057
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
  \code{\link{UCRS.5b}}, \code{\link{UCRS.5c}} where 4 indicates a lower limit fixed at 0 (not estimated), 
  5 indicates lower and upper limits being estimated (no constraints), and a, b and c indicate 
  pre-specified alpha values 1, 0.5 and 0.25, respectively.
}

\examples{

## Estimating CRS model with a lower limit equal to 0 and alpha fixed at 1
lettuce.crs.4a <- drm(weight ~ conc, data = lettuce, fct = CRS.4a())

summary(lettuce.crs.4a)

predict(lettuce.crs.4a, se.fit = TRUE)
predict(lettuce.crs.4a, interval = "confidence")
predict(lettuce.crs.4a, interval = "prediction")

plot(lettuce.crs.4a)

# Example from Sweeney et al. (2026) with U-shaped hormesis data
concVec <- c(0,0,0,1,1,1,2,2,2,4,4,4,10,10,10, 50, 50, 50)
respVec <- c(2.1, 2.0, 1.8, 1.5, 1.2, 1.1, 2.2, 1.8, 2.1, 
4.6, 4.8, 5.2, 9.8, 8.8, 10.1, 11.1, 11.3, 10.9)

ex.ucrs.5b <- drm(respVec ~ concVec, fct = UCRS.5b())
summary(ex.ucrs.5b)
# estimated parameter f: 2.62248 (0.88486), p = 0.01098

plot(ex.ucrs.5b, bp = 0.001)
fitted(ex.ucrs.5b)

predict(ex.ucrs.5b, se.fit = TRUE)
predict(ex.ucrs.5b, interval = "confidence")

}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation u-shaped}
