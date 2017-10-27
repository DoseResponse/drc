\name{lnormal}

\alias{lnormal}
\alias{LN.2}
\alias{LN.3}
\alias{LN.3u}
\alias{LN.4}

\title{Log-normal dose-response model}

\description{
  \code{lnormal} and the accompanying convenience functions provide a general framework for specifying 
  the mean function of the decreasing or incresing log-normal dose-response model.
}

\usage{
  lnormal(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), 
  method = c("1", "2", "3", "4"), ssfct = NULL,
  fctName, fctText, loge = FALSE)
  
  LN.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
  
  LN.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  LN.3u(upper = 1, fixed = c(NA, NA, NA), names = c("b", "c", "e"), ...)
  
  LN.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f 
  (see under 'Details' for the precise meaning of each parameter).}
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{character string used internally by convenience functions (optional).}  
  \item{fctText}{character string used internally by convenience functions (optional).}  
  \item{loge}{logical indicating whether or not ED50 or log(ED50) should be a parameter in the model. 
  By default ED50 is a model parameter.}
  
  \item{upper}{numeric specifying the upper horizontal asymptote in the convenience function. 
  The default is 1.}
  \item{...}{additional arguments to be passed from the convenience functions to \code{lnormal}.}
}

\details{
  For the case where log(ED50), denoted \eqn{e} in the equation below, is a parameter in the model, 
  the mean function is:
  
  \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-e)))}

  and the mean function is:
  
  \deqn{f(x) = c + (d-c)(\Phi(b(\log(x)-\log(e))))}  

  in case ED50, which is also denoted \eqn{e}, is a parameter in the model. If the former model is fitted
  any estimated ED values will need to be back-transformed subsequently in order to obtain effective doses
  on the original scale.

  The mean functions above yield the same models as those described by Bruce and Versteeg (1992),
  but using a different parameterisation (among other things the natural logarithm is used).
  
  For the case \eqn{c=0} and \eqn{d=1}, the log-normal model reduces the classic probit model (Finney, 1971) 
  with log dose as explanatory variable (mostly used for quantal data). This special case is available through
  the convenience function \code{LN.2}.
  
  The case \eqn{c=0} is available as the function \code{LN.3}, whereas the \code{LN.3u} corresponds to the special
  case where the upper horizontal asymptote is fixed (default is 1). The full four-parameter model is available
  through \code{LN.4}.
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.
}

\references{
  Finney, D. J. (1971) \emph{Probit analysis}, London: Cambridge University Press.

  Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for modeling continuous toxicity data,
  \emph{Environ. Toxicol. Chem.}, \bold{11}, 1485--1494.
}

\author{Christian Ritz}

\note{
  The function is for use with the function \code{\link{drm}}, but typically the convenience functions
  \code{link{LN.2}}, \code{link{LN.3}}, \code{link{LN.3u}}, and \code{link{LN.4}} should be used.
}

\seealso{
  The log-logistic model (\code{\link{llogistic}}) is very similar to the log-normal model at least in the middle,
  but they may differ in the tails and thus provide different estimates of low effect concentrations EC/ED.
}

%\examples{}

\keyword{models}
\keyword{nonlinear}

\concept{log-normal probit}