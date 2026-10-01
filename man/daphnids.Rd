\name{daphnids}

\alias{daphnids}

\docType{data}

\title{Daphnia test}

\description{
  The number of immobile daphnids --in contrast to mobile daphnids-- out of a total of 20 daphnids was counted 
  for several concentrations of a toxic substance.
}

\usage{data(daphnids)}

\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{dose}}{a numeric vector}
    \item{\code{no}}{a numeric vector}
    \item{\code{total}}{a numeric vector}
    \item{\code{time}}{a factor with levels \code{24h} \code{48h}}
  }
}

\details{
  The same daphnids were counted at 24h and later again at 48h.
}

\source{
  Nina Cedergreen, Faculty of Life Sciences, University of Copenhagen, Denmark.
}

%\references{}

\examples{

## Fitting a model with different parameters
## for different curves
daphnids.m1 <- drm(no/total~dose, time, weights = total, 
data = daphnids, fct = LL.2(), type = "binomial")

## Goodness-of-fit test
modelFit(daphnids.m1)

## Summary of the data
summary(daphnids.m1)

## Fitting a model with a common intercept parameter
daphnids.m2 <- drm(no/total~dose, time, weights = total, 
data = daphnids, fct = LL.2(), type = "binomial", 
pmodels = list(~1, ~time))

}
\keyword{datasets}
