\name{logistic}

\alias{logistic}

\alias{L.3}
\alias{L.4}
\alias{L.5}

%\alias{boltzmann}

%\alias{B.3}
%\alias{B.4}
%\alias{B.5}
%\alias{b3}
%\alias{b4}
%\alias{b5}

\title{The logistic model}

\description{
  The general asymmetric five-parameter logistic model for describing dose-response relationships.
}

\usage{
  logistic(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"), ssfct = NULL, 
  fctName, fctText) 

  L.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  L.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
  L.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)  

%  boltzmann(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), fctName, fctText)
%
%  B.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"))
%  B.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"))
%  B.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"))  
}


\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The order of the parameters is: b, c, d, e, f (see under 'Details').}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}    
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}
  \item{...}{Additional arguments (see \code{\link{llogistic}}).}    
}

\details{
  The default arguments yields the five-parameter logistic mean function given by the expression
  
  \deqn{ f(x) = c + \frac{d-c}{(1+\exp(b(x - e)))^f}}
  
  The model is different from the log-logistic models \code{\link{llogistic}} and \code{\link{llogistic2}}
  where the term \deqn{log(x)} is used instead of \deqn{x}.
  
  The model is sometimes referred to as the Boltzmann model.
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}

% \references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

## Fitting the four-parameter logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.4())
summary(ryegrass.m1)

## Fitting an asymmetric logistic model
##  requires installing the package 'NISTnls'
# Ratkowsky3.m1 <- drm(y~x, data = Ratkowsky3, 
# fct = L.5(fixed = c(NA, 0, NA, NA, NA)))
# plot(Ratkowsky3.m1)
# summary(Ratkowsky3.m1)  
## okay agreement with NIST values
##  for the two parameters that are the same

}
\keyword{models}
\keyword{nonlinear}
