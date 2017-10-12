\name{fitted.drc}

\alias{fitted.drc}

\title{Extract fitted values from model}

\description{
  Extracts fitted values from an object of class 'drc'. 
}

\usage{

  \method{fitted}{drc}(object, ...)

}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
}

\value{
  Fitted values extracted from 'object'.
}

\author{Christian Ritz}

\examples{

ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
plot(fitted(ryegrass.m1), residuals(ryegrass.m1))  # a residual plot

}
\keyword{models}
\keyword{nonlinear}
