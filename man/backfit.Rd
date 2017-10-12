\name{backfit}

\alias{backfit}

\title{
  Calculation of backfit values from a fitted dose-response model
}

\description{
  By inverse regression backfitted dose values are calculated for the mean response per dose.
}

\usage{
backfit(drcObject)
}

\arguments{
  \item{drcObject}{an object of class 'drc'.}
}                      

\value{
  Two columns with the original dose values and the corresponding backfitted values using the fitted dose-response model. 
  For extreme dose values (e.g., high dose ) the backfitted values may not be well-defined (see the example below). 
}

\references{
  ??
}

\author{
  Christian Ritz after a suggestion from Keld Sorensen.
}

\seealso{
  A related function is \code{\link{ED.drc}}.
}

\examples{

ryegrass.LL.4 <- drm(rootl~conc, data=ryegrass, fct=LL.4())

backfit(ryegrass.LL.4)

}

\keyword{models}
\keyword{nonlinear}
