\name{NEC}

\Rdversion{1.1}

\alias{NEC}
\alias{NEC.2}
\alias{NEC.3}
\alias{NEC.4}

\title{
  Dose-response model for estimation of no effect concentration (NEC).
}

\description{
  The no effect concentration has been proposed as an alternative to both the classical no observed effect concentration
  (NOEC) and the regression-based EC/ED approach. The NEC model is a dose-response model with a threshold below 
  which the response is assumed constant and equal to the control response.
}

\usage{
  NEC(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), fctName, fctText)

  NEC.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
  
  NEC.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
  
  NEC.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)

}

\arguments{
  \item{fixed}{numeric vector specifying which parameters are fixed and at what value they are fixed. 
  NAs are used for parameters that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").
  The default is reasonable (see under 'Usage').}
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}
  \item{upper}{numeric value. The fixed, upper limit in the model. Default is 1.}
  \item{...}{additional arguments in \code{\link{NEC}}}  
}


\details{
  The NEC model function proposed by Pires \emph{et al} (2002) is defined as follows
  
  \deqn{ f(x) = c + (d-c)\exp(-b(x-e)I(x-e)) + \frac{d2}{1+\exp(b2(\log(x)-\log(e2)))}}
  
  where \eqn{I(x-e)} is the indicator function. It is equal to 0 for \eqn{x<=e} and equal 1 for \eqn{x>e}.
  
  In other words: The parameter e in \code{NEC} in "drc" corresponds to the parameter c' in Pires \emph{et al} (2002),
  the parameter b in \code{NEC} in "drc" corresponds to the parameter m' in Pires \emph{et al} (2002), the parameter d 
  in \code{NEC} in "drc" corresponds to the parameter l' in Pires \emph{et al} (2002), and finally the parameter c in 
  \code{NEC} in "drc" (the lower horizontal limit) is (implictly) fixed at 0 in Pires \emph{et al} (2002)
}

\value{
  The value returned is a list containing the nonlinear function, the self starter function
  and the parameter names.
}
\references{
  
  Pires, A. M., Branco, J. A., Picado, A., Mendonca, E. (2002)
  Models for the estimation of a 'no effect concentration',
  \emph{Environmetrics}, \bold{13}, 15--27.

}
\author{
  Christian Ritz
}

%\note{}

%\seealso{}

\examples{

nec.m1 <- drm(rootl~conc, data=ryegrass, fct=NEC.4())

summary(nec.m1)

plot(nec.m1)

abline(v=coef(nec.m1)[4], lty=2)  # showing the estimated threshold 

}

\keyword{models}
\keyword{nonlinear}
