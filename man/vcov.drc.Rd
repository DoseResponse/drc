\name{vcov.drc}

\alias{vcov.drc}

\title{Calculating variance-covariance matrix for objects of class 'drc'}

\description{
  'vcov' returns the estimated variance-covariance matrix for the parameters in the non-linear function.
}

\usage{
  \method{vcov}{drc}(object, ..., corr = FALSE, od = FALSE, pool = TRUE, unscaled = FALSE)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{corr}{logical. If TRUE a correlation matrix is returned.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used. This argument only makes a difference for
  binomial data.}
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.} 
  \item{unscaled}{logical. If TRUE the unscaled variance-covariance is returned. This argument only makes a difference
  for continuous data.}
}

\value{
  A matrix of estimated variances and covariances.
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
vcov(ryegrass.m1)
vcov(ryegrass.m1, corr = TRUE)

}
\keyword{models}
\keyword{nonlinear}
