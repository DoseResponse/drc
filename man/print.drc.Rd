\name{print.drc}

\alias{print.drc}

\title{Printing key features}

\description{
  'print' displays brief information on an object of class 'drc'.
}

\usage{

  \method{print}{drc}(x, ..., digits = max(3, getOption("digits") - 3))

}

\arguments{
  \item{x}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{digits}{an integer giving the number of digits of the parameter coefficients. Default is 3.} 
}

\author{Christian Ritz}

\examples{

## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())

## Displaying the model fit
print(ryegrass.m1)
ryegrass.m1  # gives the same output as the previous line

}
\keyword{models}
\keyword{nonlinear}
