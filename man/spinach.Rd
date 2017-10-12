\name{spinach}

\alias{spinach}

\docType{data}

\title{Inhibition of photosynthesis}

\description{
  Data from an experiment investigating the inhibition of photosynthesis in response to two synthetic 
  photosystem II inhibitors, the herbicides diuron and bentazon.
  More specifically, the effect of oxygen consumption of thylakoid membranes (chloroplasts) from spinach
  was measured after incubation with the synthetic inhibitors in 5 assays, resulting in 5 dose-response curves.
}

\usage{data(spinach)}

\format{
  A data frame with 105 observations on the following four variables:
  \describe{
    \item{CURVE}{a numeric vector specifying the assay or curve (a total of 5 independent assays where used in this experiment).}
    \item{HERBICIDE}{a character vector specifying the herbicide applied: bentazon or diuron.}
    \item{DOSE}{a numeric vector giving the herbicide concentration in muMol.}
    \item{SLOPE}{a numeric vector with the measured response: oxygen consumption of thylakoid membranes.}
  }
}

\details{
  The experiment is described in more details by Streibig (1998).
}

\source{
  Streibig, J. C. (1998) Joint action of natural and synthetic photosystem II inhibitors, \emph{Pesticide Science}, \bold{55}, 137--146.
}

\examples{

## Displaying the first rows in the dataset
head(spinach)  # displaying first 6 rows in the data set

}
\keyword{datasets}
