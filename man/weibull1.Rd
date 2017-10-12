\name{weibull1}

\alias{weibull1}
\alias{weibull2}
\alias{weibull2x}

\title{Weibull model functions}

\description{
  'weibull' and 'weibull2' provide a very general way of specifying Weibull dose response functions,
  under various constraints on the parameters.
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
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}                    
}

\details{
  As pointed out in Seber and Wild (1989), there exist two different parameterisations of the Weibull model. They
  do not yield the same fitted curve for a given dataset (see under Examples).

  The four-parameter Weibull type 1 model ('weibull1') is
  \deqn{ f(x) = c + (d-c) \exp(-\exp(b(\log(x)-\log(e)))).}
 
  Thw four-parameter Weibull type 2 model ('weibull2') is
  \deqn{ f(x) = c + (d-c) (1 - \exp(-\exp(b(\log(x)-\log(e))))).}  
  
  Both four-parameter model functions are asymmetric with inflection point at the dose equal \eqn{e}. 
}

\value{
  The value returned is a list containing the non-linear function, the self starter function
  and the parameter names.  
}

\references{
  Seber, G. A. F. and Wild, C. J (1989) 
  \emph{Nonlinear Regression}, 
  New York: Wiley \& Sons (pp. 338--339).
}

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{drm}}.

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
