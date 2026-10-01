\name{heartrate}

\alias{heartrate}

\docType{data}

\title{Heart rate baroreflexes for rabbits}

\description{
  The dataset contains measurements of mean arterial pressure (mmHG) and heart rate (b/min) for a baroreflex curve.
}

\usage{data(heartrate)}

\format{
  A data frame with 18 observations on the following 2 variables.
  \describe{
    \item{\code{pressure}}{a numeric vector containing measurements of arterial pressure.}
    \item{\code{rate}}{a numeric vector containing measurements of heart rate.}
  }
}

\details{
  The dataset is an example of an asymmetric dose-response curve, that is not
  easily handled using the log-logistic or Weibull models (\code{\link{LL.4}}, \code{\link{LL.5}}, 
  \code{\link{W1.4}} and \code{\link{W2.4}}), whereas the \code{baro5} model provides a nice fit.
}

\source{
  Ricketts, J. H. and Head, G. A. (1999) A five-parameter logistic equation for investigating asymmetry of
  curvature in baroreflex studies, 
  \emph{Am. J. Physiol. (Regulatory Integrative Comp. Physiol. 46)}, \bold{277}, 441--454.
}

\examples{

## Fitting the baro5 model
heartrate.m1 <- drm(rate~pressure, data=heartrate, fct=baro5())
plot(heartrate.m1)

coef(heartrate.m1)

#Output:
#b1:(Intercept) b2:(Intercept)  c:(Intercept)  d:(Intercept)  e:(Intercept)
#      11.07984       46.67492      150.33588      351.29613       75.59392

## Inserting the estimated baro5 model function in deriv()
baro5Derivative <- deriv(~ 150.33588 + ((351.29613 - 150.33588)/
(1 + (1/(1 + exp((2 * 11.07984 * 46.67492/(11.07984 + 46.67492)) * 
(log(x) - log(75.59392 ))))) * (exp(11.07984 * (log(x) - log(75.59392)))) + 
(1 - (1/(1 + exp((2 * 11.07984 * 46.67492/(11.07984 + 46.67492)) * 
(log(x) - log(75.59392 )))))) * (exp(46.67492 * (log(x) - log(75.59392 )))))), "x", function(x){})

## Plotting the derivative
#pressureVector <- 50:100
pressureVector <- seq(50, 100, length.out=300)
derivativeVector <- attr(baro5Derivative(pressureVector), "gradient")
plot(pressureVector, derivativeVector, type = "l")

## Finding the minimum
pressureVector[which.min(derivativeVector)]


}
\keyword{datasets}
