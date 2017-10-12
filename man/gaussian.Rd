\name{gaussian}

\alias{gaussian}
\alias{lgaussian}


\title{
  Normal and log-normal biphasic dose-response models
}

\description{
  Model functions for fitting symmetric or skewed bell-shaped/biphasic dose-response patterns.
}

\usage{
    gaussian(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
    method = c("1", "2", "3", "4"), ssfct = NULL, fctName, fctText, loge = FALSE)
    
    lgaussian(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", 
    "d", "e", "f"), method = c("1", "2", "3", "4"), ssfct = NULL, 
    fctName, fctText, loge = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
    NAs for parameter that are not fixed.}
    
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f 
    (see under 'Details').}
  
  \item{method}{character string indicating the self starter function to use.} 
  
  \item{ssfct}{a self starter function to be used.}  
  
  \item{fctName}{optional character string used internally by convenience functions.}  
  
  \item{fctText}{optional character string used internally by convenience functions.} 
  
  \item{loge}{logical indicating whether or not e or log(e) should be a parameter in the model. By default e is a model parameter.}

}

\details{
  Details yet to be provided.
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{drm}}.
}

\keyword{models}
\keyword{nonlinear}