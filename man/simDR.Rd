\name{simDR}

\alias{simDR}

\title{Simulating ED values under various scenarios}

\description{
  Simulating ED values for a given model and given dose values.
}

\usage{
  simDR(mpar, sigma, fct, noSim = 1000, conc, edVec = c(10, 50), seedVal = 20070723)
}

\arguments{
  \item{mpar}{numeric vector of model parameters}
  \item{sigma}{numeric specifying the residual standard deviation}
  \item{fct}{list supplying the chosen mean function}
  \item{conc}{numeric vector of concentration/dose values}
  \item{edVec}{numeric vector of ED values to estimate in each simulation}
  \item{noSim}{numeric giving the number of simulations}
  \item{seedVal}{numeric giving the seed used to initiate the random number generator}
}

\details{
  The arguments \code{mpar} and \code{sigma} are typically obtained from a previous model fit.

  Only dose-response models assuming normally distributed errors can be used.
}

\value{
  A list of matrices with as many components as there are chosen ED values. The entries in the matrices are
  empirical standard deviations of the estimated ED values. Row-wise from top to bottom more and more
  concentration/dose values are included in the simulations; top row starting with 5 concentrations. The
  number of replicates increases column by column from left to right. 
  
  The list is returned invisbly as the matrices also are displayed.
}

%\references{ ~put references to the literature/web site here ~ }

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

ryegrass.m1 <- drm(ryegrass, fct=LL.4())

simDR(coef(ryegrass.m1), sqrt(summary(ryegrass.m1)$resVar), LL.4(), 2, 
c(1.88, 3.75, 7.50, 0.94, 15, 0.47, 30, 0.23, 60), seedVal = 200710291)

}
\keyword{models}
\keyword{nonlinear}
