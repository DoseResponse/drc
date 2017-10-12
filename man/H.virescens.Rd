\name{H.virescens}

\alias{H.virescens}

\docType{data}

\title{Mortality of tobacco budworms}

\description{
  For three days, moths of the tobacco budworm (\emph{Heliothis virescens}) were exposed 
  to doses of the pyrethroid trans-cypermethrin.
}

\usage{data(H.virescens)}

\format{
  A data frame with 12 observations on the following 4 variables.
  \describe{
    \item{\code{dose}}{a numeric vector of dose values (\eqn{\mu g})}
    \item{\code{numdead}}{a numeric vector of dead or knocked-down moths}
    \item{\code{total}}{a numeric vector of total number of moths}
    \item{\code{sex}}{a factor with levels \code{F} \code{M} denoting a grouping according to sex}
  }
}

\details{
  In Venables and Ripley (2002), these data are analysed using a logistic regression with base-2 logarithm of dose 
  as explanatory variable.
}

\source{
  Venables, W. N. and Ripley, B. D (2002) \emph{Modern Applied Statistics with S}, New York: Springer (fourth edition).
}

%\references{}

\examples{

## Fitting dose-response model (log-logistic with common slope)
Hv.m1 <- drm(numdead/total~dose, sex, weights = total, data = H.virescens, fct = LL.2(), 
pmodels = list(~ 1, ~ sex - 1), type = "binomial")
summary(Hv.m1)

## Fitting the same model as in Venables and Riply (2002)
Hv.m2 <- glm(cbind(numdead, total-numdead) ~ sex + I(log2(dose)) - 1, data = H.virescens, 
family = binomial)

## Comparing the fits
logLik(Hv.m1)
logLik(Hv.m2)

## Estimated ED values (matching those given in MASS)
ED(Hv.m1, c(25, 50, 75))


}
\keyword{datasets}
