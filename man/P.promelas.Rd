\name{P.promelas}

\alias{P.promelas}

\docType{data}

\title{Effect of sodium pentachlorophenate on growth of fathead minnow}

\description{
  Fathead minnows (\emph{Pimephales promelas}) were exposed to sodium pentachlorophenate
  concentrations ranging from 32 to 512 micro g/L in a 7-day larval survival and growth test.
  The average dry weight was measured.
}

\usage{data(P.promelas)}

\format{
  A data frame with 24 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of sodium pentachlorophenate concentrations (micro g/L).}
    \item{\code{dryweight}}{a numeric vector dry weights (mg)}
  }
}

\details{
  The data are analysed in Bruce and Versteeg (1992) using a log-normal 
  dose-response model (using the logarithm with base 10).
}

\source{
  Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for modeling continuous toxicity data,
  \emph{Environ. Toxicol. Chem.}, \bold{11}, 1485--1494.
}

%\references{}

\examples{

## Model with ED50 on log scale as parameter
p.prom.m1<-drm(dryweight~conc, data=P.promelas, fct=LN.3())

plot(fitted(p.prom.m1), residuals(p.prom.m1))

plot(p.prom.m1, type="all", broken=TRUE, xlim=c(0,1000))
summary(p.prom.m1)
ED(p.prom.m1, c(10,20,50), interval="delta")

## Model with ED50 as parameter
p.prom.m2<-drm(dryweight~conc, data=P.promelas, fct=LN.3(loge=TRUE))
summary(p.prom.m2)
ED(p.prom.m2, c(10,20,50), interval="fls")

}
\keyword{datasets}
