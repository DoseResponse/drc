\name{searchdrc}
\alias{searchdrc}

\title{Searching through a range of initial parameter values to obtain convergence}

\description{
  'searchdrc' provides a facility for searching through a range of parameter values (one-dimensional)
  in order to obtain convergence of the estimation procedure.
}

\usage{
searchdrc(object, which, range, len = 50)
}

\arguments{
  \item{object}{an object of class 'drc'. The object can be from a model that could not fitted.}
  \item{which}{a character string containing the parameter name}
  \item{range}{a numeric vector of length 2 specifying the interval endpoints for the range.}
  \item{len}{numeric. The number of points in the interval.}
}

\details{
  The function goes through the range with increments such that in total at most \code{len} sets of parameter values
  are used as initial values for the estimation procedure. You would need to identify the parameter which is most likely to
  cause problems for the estimation procedure.
}

\value{
  An object of class 'drc'. 
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{ ~~further notes~~ }

%\seealso{ ~~objects to See Also as \code{\link{~~fun~~}}, ~~~ }

\examples{

## No example yet

}
\keyword{models}
\keyword{nonlinear}
