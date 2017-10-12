\name{twophase}

\Rdversion{1.1}

\alias{twophase}

\title{
  Two-phase dose-response model
}

\description{
   The two-phase dose-response model is a combination of log-logistic models that should be useful for describing
   more complex dose-response patterns.
}

\usage{
  twophase(fixed = c(NA, NA, NA, NA, NA, NA, NA), 
  names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2"), fctName, fctText)
}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at what value they are fixed. 
  NAs are used for parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage').}
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}
}

\details{
  Following Groot \emph{et al} (1996) the two-phase model function is defined as follows
  
  \deqn{ f(x) = c + \frac{d1-c}{1+\exp(b1(\log(x)-\log(e1)))} + \frac{d2}{1+\exp(b2(\log(x)-\log(e2)))}}
  
  For each of the two phases, the parameters have the same interpretation as in the ordinary log-logistic
  model.
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}

\references{
  Groot, J. C. J., Cone, J. W., Williams, B. A., Debersaques, F. M. A., Lantinga, E. A. (1996)
  Multiphasic analysis of gas production kinetics for in vitro fermentation of ruminant feeds,
  \emph{Animal Feed Science Technology}, \bold{64}, 77--89.  
}

\author{
  Christian Ritz
}

%\note{}

\seealso{
  The basic component in the two-phase model is the log-logistic model
  \code{\link{llogistic}}.
}

%\examples{}

\keyword{models}
\keyword{nonlinear}
