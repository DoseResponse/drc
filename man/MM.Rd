\name{MM}

\alias{MM.2}
\alias{MM.3}

\title{Michaelis-Menten model}

\description{
  The functions can be used to fit (shifted) Michaelis-Menten models that are used
  for modeling enzyme kinetics, weed densities etc.
}

\usage{ 
  MM.2(fixed = c(NA, NA), names = c("d", "e"), ...)
  
  MM.3(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)  
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":").}
  \item{...}{additional arguments from convenience functions to \code{\link{llogistic}}.}
}

\details{
  The model is defined by the three-parameter model function
  
  \deqn{f(x, (c, d, e)) = c + \frac{d-c}{1+(e/x)}}
  
  It is an increasing as a function of the dose \eqn{x}, attaining the lower limit \eqn{c} at dose 0 (\eqn{x=0}) 
  and the upper limit \eqn{d} for infinitely large doses. The parameter \eqn{e} corresponds to the dose yielding a response
  halfway between \eqn{c} and \eqn{d}. 
  
  The common two-parameter Michaelis-Menten model (\code{MM.2}) is obtained by
  setting \eqn{c} equal to 0.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{}

\author{Christian Ritz}

\note{
  At the moment the implementation cannot deal with infinite concentrations.
}

\seealso{Related models are the asymptotic regression models \code{\link{AR.2}} and \code{\link{AR.3}}.}

\examples{

## Fitting Michaelis-Menten model
met.mm.m1 <- drm(gain~dose, product, data=methionine, fct=MM.3(), 
pmodels = list(~1, ~factor(product), ~factor(product)))
plot(met.mm.m1, log = "", ylim=c(1450, 1800))
summary(met.mm.m1)
ED(met.mm.m1, c(10, 50))

## Calculating bioefficacy: approach 1
coef(met.mm.m1)[4] / coef(met.mm.m1)[5] * 100

## Calculating bioefficacy: approach 2
EDcomp(met.mm.m1, c(50,50))

## Simplified models
met.mm.m2a <- drm(gain~dose, product, data=methionine, fct=MM.3(), 
pmodels = list(~1, ~factor(product), ~1))
anova(met.mm.m2a, met.mm.m1)  # model reduction not possible

met.mm.m2b <- drm(gain~dose, product, data=methionine, fct=MM.3(), 
pmodels = list(~1, ~1, ~factor(product)))
anova(met.mm.m2b, met.mm.m1)  # model reduction not possible

}
\keyword{models}
\keyword{nonlinear}
