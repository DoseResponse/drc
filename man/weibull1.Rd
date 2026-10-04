\name{weibull1}

\alias{weibull1}
\alias{weibull2}
\alias{weibull2x}

\title{Weibull model functions}

\description{
  'weibull1' and 'weibull2' provide a very general ways of specifying Weibull type 1 and 2 dose-response models,
  allowing for various constraints on the model parameters.
}

\usage{
  weibull1(fixed = c(NA, NA, NA, NA), 
           names = c("b", "c", "d", "e"),
           method = c("1", "2", "3", "4"), 
           ssfct = NULL,
           fctName, fctText)
 
  weibull2(fixed = c(NA, NA, NA, NA), 
           names = c("b", "c", "d", "e"),
           method = c("1", "2", "3", "4"), 
           ssfct = NULL,
           fctName, fctText)

  weibull2x(fixed = rep(NA, 5), 
           names = c("b", "c", "d", "e", "t0"),
           method = c("1", "2", "3", "4"), 
           ssfct = NULL,
           fctName, fctText)

}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e (see under 'Details').}
  \item{method}{character string indicating the self starter function to use; default is "1" (see under Details).}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}                    
}

\details{
  As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They
  do not yield the same fitted curve for a given dataset (see under Examples).

  The four-parameter Weibull type 1 model ('weibull1') is defined by the model function:
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
 
  Likewise, the four-parameter Weibull type 2 model ('weibull2') is defined by the model function:
  \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).}  
  
  Both model functions result in asymmetric dose-response curves. Weibull type 2 exhibits 
  a rapid descent from the upper asymptote and a slow approach to the lower asymptote (in case of a decreasing response). 
  Weibull type 1 has the opposite behaviour.
  
  The four-parameter Weibull type 1 and Weibull type 2 model functions have an inflection point at the dose equal 
  to \eqn{e(1-1/b)^{1/b}}, which is always smaller, but often not much smaller, than the model parameter \eqn{e}, 
  which is not the EC/ED50. 

  The default self starter function (method = "1") is based on linear regression applied to the linearized response obtained
  using a suitable log-log transformation. Methods "2" and "4" are refinements of method "1" (Normolle, 1993). Method "3" 
  is based on increments in the response.
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.  
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) 
  \emph{Nonlinear Regression}, 
  New York: Wiley and Sons (pp. 338--339).

  Normolle, D.P. (1993) 
  An algorithm for robust non-linear analysis of radioimmunoassays and other bioassays,
  \emph{Statistics in Medicine}, \bold{12}, 2025--2042. 
  https://doi.org/10.1002/sim.4780122106

  Ritz, C. and Jensen, S. M. and Gerhard, D. and Streibig, J. C. (2019)
  \emph{Dose-Response Analysis using R}, 
  Boca Raton: CRC Press (pp. 184--186).
}

\author{Christian Ritz}

\note{
  These model functions are for use with the function \code{\link{drm}}.
}

\seealso{
  For convenience several special cases of the function 'weibull1' are available: 
  \code{\link{W1.2}}, \code{\link{W1.3}} and \code{\link{W1.4}}.  

  Special cases of 'weibull2' are: 
  \code{\link{W2.2}}, \code{\link{W2.3}} and \code{\link{W2.4}}.  
  
  These convenience functions should be used rather than the underlying functions 
  \code{weibull1} and \code{weibull2}.
}

\examples{

## Fitting two different Weibull models
ryegrass.m1 <- drm(ryegrass, fct = W1.4())
plot(ryegrass.m1, conLevel=0.5)

ryegrass.m2 <- drm(ryegrass, fct = W2.4())
plot(ryegrass.m2, conLevel=0.5, add = TRUE, type = "none", col = 2)
# you could also look at the ED values to see the difference

## A four-parameter Weibull model with b fixed at 1
ryegrass.m3 <- drm(ryegrass, fct = W1.4(fixed = c(1, NA, NA, NA)))
summary(ryegrass.m3)

## A four-parameter Weibull model with the constraint b>3
ryegrass.m4 <- drm(ryegrass, fct = W1.4(), lowerl = c(3, -Inf, -Inf, -Inf), 
control = drmc(constr=TRUE))
summary(ryegrass.m4)

}
\keyword{models}
\keyword{nonlinear}
