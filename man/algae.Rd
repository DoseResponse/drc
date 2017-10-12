\name{algae}

\alias{algae}

\docType{data}

\title{Volume of algae as function of increasing concentrations of a herbicide}

\description{
  Dataset from an experiment exploring the effect of increasing concentrations of a herbicide on 
  the volume of the treated algae.
}

\usage{data(algae)}

\format{
  A data frame with 14 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of concentrations.}
    \item{\code{vol}}{a numeric vector of response values, that is relative change in volume.}
  }
}

\details{
  This datasets requires a cubic root transformation in order to stabilise the variance.
}

\source{
  Meister, R. and van den Brink, P. (2000)
  \emph{The Analysis of Laboratory Toxicity Experiments},  
  Chapter 4 in \emph{Statistics in Ecotoxicology}, Editor: T. Sparks, 
  New York: John Wiley \& Sons, (pp. 114--116). 
}

%\references{}

\examples{

algae.m1 <- drm(vol~conc, data=algae, fct=LL.3())
summary(algae.m1)

algae.m2 <- boxcox(algae.m1)
summary(algae.m2)

}
\keyword{datasets}
