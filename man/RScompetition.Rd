\name{RScompetition}

\alias{RScompetition}

\docType{data}

\title{Competition between two biotypes}

\description{
  To assess the competitive ability between two biotypes of \emph{Lolium rigidum}, one resistant to glyphosate
  and the other a sensitive wild type, the density of resistant and sensitive biotypes was counted after
  germination.
}

\usage{data(RScompetition)}

\format{
  A data frame with 49 observations on the following 3 variables.
  \describe{
    \item{\code{z}}{a numeric vector with densities of the resistant biotype (plants/m2)}
    \item{\code{x}}{a numeric vector with densities of the sensitive biotype (plants/m2)}
    \item{\code{biomass}}{a numeric vector of biomass weight (g/plant)}
  }
}

\details{
  A hyperbolic model (Jensen, 1993) is describing the data reasonably well.
}

\source{
  The dataset is from Pedersen et al (2007).
}

\references{
  Jensen, J. E. (1993) Fitness of herbicide-resistant weed biotypes described by competition models, 
  \emph{Proceedings of the 8th EWRS Symposium, 14-16 June, Braunschweig, Germany}, 
  \bold{1}, 25--32.

  Pedersen, B. P. and Neve, P. and Andreasen, C. and Powles, S. (2007) Ecological fitness of a glyphosate 
  resistant \emph{Lolium rigidum} population: Growth and seed production along a competition gradient,
  \emph{Basic and Applied Ecology}, \bold{8}, 258--268.
}

%\examples{}

\keyword{datasets}
