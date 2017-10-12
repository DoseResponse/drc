\name{nasturtium}

\alias{nasturtium}

\docType{data}

\title{Dose-response profile of degradation of agrochemical using nasturtium}

\description{
  Estimation of the degradation profile of an agrochemical based on soil samples at depth 0-10cm
  from a calibration experiment.
}

\usage{data(nasturtium)}

\format{
  A data frame with 42 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of concentrations (g/ha)}
    \item{\code{weight}}{a numeric vector of plant weight (mg) after 3 weeks' growth}
  }
}

\details{
  It is an experiment with seven concentrations and six replicates per concentration. \emph{Nasturtium}
  is sensitive and its weight reduces noticeable at low concentrations.

  Racine-Poon (1988) suggests using a three-parameter log-logistic model.
}

\source{
  Racine-Poon, A. (1988) A Bayesian Approach to Nonlinear Calibration Problems, 
  \emph{J. Am. Statist. Ass.}, \bold{83}, 650--656.
}

%\references{}

\examples{

nasturtium.m1 <- drm(weight~conc, data=nasturtium, fct = LL.3())

modelFit(nasturtium.m1)

plot(nasturtium.m1, type = "all", log = "", xlab = "Concentration (g/ha)", ylab = "Weight (mg)")
}

\keyword{datasets}
