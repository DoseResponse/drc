\name{logLik.drc}

\alias{logLik.drc}

\title{Extracting the log likelihood}

\description{
  \code{loglik} extracts the value of the log likelihood function evaluated at the parameter estimates.
}
\usage{
  \method{logLik}{drc}(object, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
   The evaluated log likelihood as a numeric value and the corresponding degrees of freedom as well as the number of observations as attributes.
}

\note{
    The value of the log likelihood could be used to compare model fits of the same data based on different 
    dose-response models or based on the same model but fitted different algorithms, software programmes, or 
    starting values. For comparisons: Larger is better.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4()) 
logLik(ryegrass.m1)

}
\keyword{models}
\keyword{nonlinear}
