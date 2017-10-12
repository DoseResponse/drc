\name{PR}

\alias{PR}

\title{Expected or predicted response}

\description{
  The function returns the expected or predicted response for specified dose values.
}

\usage{
  PR(object, xVec, ...)
}

\arguments{
  \item{object}{object of class \code{drc} obtaining fitting a dose-response model.}
  \item{xVec}{numeric vector of dose values.}
  \item{\dots}{additional arguments to be supplied to \code{\link[drc]{predict.drc}}. No effect at the moment.}
}

\details{
  This function is a convenience function for easy access to predicted values.
}

\value{
  A numeric vector of predicted values or possibly a matrix of predicted values and corresponding standard errors.
}

%\references{}

\author{
  Christian Ritz after a suggestion from Andrew Kniss.
}

%\note{ ~~further notes~~ }

\seealso{Predictions can also be obtained using \code{\link[drc]{predict.drc}}.}

\examples{

ryegrass.m1 <- drm(ryegrass, fct = LL.4())
PR(ryegrass.m1, c(5, 10))

ryegrass.m2 <- drm(ryegrass, fct = LL2.4())
PR(ryegrass.m2, c(5, 10))

spinach.m1 <- drm(SLOPE~DOSE, CURVE, data=spinach, fct = LL.4())
PR(spinach.m1, c(5, 10))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{nonlinear}
