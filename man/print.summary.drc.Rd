\name{print.summary.drc}

\alias{print.summary.drc}

\title{Printing summary of non-linear model fits}

\description{
  This method produces formatted output of the summary statistics: parameter estimates, estimated standard errors, 
  z-test statistics and corresponding p-values.  
}

\usage{

  \method{print}{summary.drc}(x, ...)

}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
  The object (argument \code{x}) is returned invisibly.
}

\author{Christian Ritz}

\examples{

ryegrass.m1 <- drm(rootl~conc, data=ryegrass, fct= LL.4())

summary(ryegrass.m1)  

}
\keyword{models}
\keyword{nonlinear}
