\name{lepidium}

\alias{lepidium}

\docType{data}

\title{Dose-response profile of degradation of agrochemical using lepidium}

\description{
  Estimation of the degradation profile of an agrochemical based on soil samples at depth 0-10cm
  from a calibration experiment.
}

\usage{data(lepidium)}

\format{
  A data frame with 42 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of concentrations (g/ha)}
    \item{\code{weight}}{a numeric vector of plant weight (g) after 3 weeks' growth}
  }
}

\details{
  It is an experiment with seven concentrations and six replicates per concentration. \emph{Lepidium}
  is rather robust as it only responds to high concentrations.
}

\source{
  Racine-Poon, A. (1988) A Bayesian Approach to Nonlinear Calibration Problems, 
  \emph{J. Am. Statist. Ass.}, \bold{83}, 650--656.
}

%\references{}

\examples{

lepidium.m1 <- drm(weight~conc, data=lepidium, fct = LL.4())

modelFit(lepidium.m1)

plot(lepidium.m1, type = "all", log = "")
}
\keyword{datasets}
