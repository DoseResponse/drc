\name{AR}

\alias{AR.2}
\alias{AR.3}

\title{Asymptotic regression model}

\description{
  Providing the mean function and the corresponding self starter function for the asymptotic regression model.
}

\usage{
  AR.2(fixed = c(NA, NA), names = c("d", "e"), ...)

  AR.3(fixed = c(NA, NA, NA), names = c("c", "d", "e"), ...)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{vector of character strings giving the names of the parameters (should not contain ":").}
  \item{...}{additional arguments to be passed from the convenience functions.}     
}

\details{
  The asymptotic regression model is a three-parameter model with mean function:
  
  \deqn{ f(x) = c + (d-c)(1-\exp(-x/e))}
  
  The parameter \eqn{c} is the lower limit (at \eqn{x=0}), the parameter \eqn{d} is the upper limit 
  and the parameter \eqn{e>0} is determining the steepness of the increase as \eqn{x}.
}

\value{
  A list of class \code{drcMean}, containing the mean function, the self starter function,
  the parameter names and other components such as derivatives and a function for calculating ED values.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

\note{
  The functions are for use with the function \code{\link{drm}}.
}

\seealso{
  A very similar, but monotonously decreasing model is the exponential decay model:
  \code{\link{EXD.2}} and \code{\link{EXD.3}}.
}

\examples{

## First model
met.as.m1<-drm(gain ~ dose, product, data = methionine, fct = AR.3(), 
pmodels = list(~1, ~factor(product), ~factor(product)))
plot(met.as.m1, log = "", ylim = c(1450, 1800))
summary(met.as.m1)

## Calculating bioefficacy: approach 1
coef(met.as.m1)[5] / coef(met.as.m1)[4] * 100

## Calculating bioefficacy: approach 2
EDcomp(met.as.m1, c(50,50))

## Simplified models
met.as.m2<-drm(gain ~ dose, product, data = methionine, fct = AR.3(), 
pmodels = list(~1, ~1, ~factor(product)))
anova(met.as.m2, met.as.m1)  # simplification not possible

met.as.m3 <- drm(gain ~ dose, product, data = methionine, fct = AR.3(), 
pmodels = list(~1, ~factor(product), ~1))
anova(met.as.m3, met.as.m1)  # simplification not possible

}
\keyword{models}
\keyword{nonlinear}

