\name{coef.drc}

\alias{coef.drc}

\title{Extract Model Coefficients}

\description{
   Extract parameter estimates.
}

\usage{

  \method{coef}{drc}(object, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
  A vector of parameter coefficients which are extracted from the model object 'object'.
}

\author{Christian Ritz}

\note{This function may work even in cases where 'summary' does not, because the parameter estimates 
      are retrieved directly from the model fit object without any additional computations of summary statistics and standard errors.}


\seealso{A more comprehensive summary is obtained using \code{\link{summary.drc}}.}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
coef(ryegrass.m1)

}
\keyword{models}
\keyword{nonlinear}
