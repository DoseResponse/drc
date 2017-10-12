\name{summary.drc}

\alias{summary.drc}

\title{Summarising non-linear model fits}

\description{
  'summary' compiles a comprehensive summary for objects of class 'drc'. 
}

\usage{

  \method{summary}{drc}(object, od = FALSE, pool = TRUE, ...)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{od}{logical. If TRUE adjustment for over-dispersion is used.}  
  \item{pool}{logical. If TRUE curves are pooled. Otherwise they are not. This argument only works for models with
  independently fitted curves as specified in \code{\link{drm}}.} 
  \item{...}{additional arguments.}
}

\value{
  A list of summary statistics that includes parameter estimates and estimated standard errors.
}

\author{Christian Ritz}

\note{
  Examples on usage are for instance found in the help pages of \code{\link{ryegrass}} and \code{\link{secalonic}}.
}

%\examples{}

\keyword{models}
\keyword{nonlinear}
