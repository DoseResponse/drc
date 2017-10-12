\name{CIcompX}
\alias{CIcompX}
\alias{CIcomp}
\alias{plotFACI}

\title{
Calculation of combination index for binary mixtures
}

\description{
For single mixture data combination indices for effective doses as well as effects may be calculated and visualized.
}

\usage{
CIcomp(mixProp, modelList, EDvec)

CIcompX(mixProp, modelList, EDvec, EDonly = FALSE)

plotFACI(effList, indAxis = c("ED", "EF"), caRef = TRUE, 
showPoints = FALSE, add = FALSE, ylim, ...)

}

\arguments{

  \item{mixProp}{a numeric value between 0 and 1 specifying the mixture proportion/ratio for
  the single mixture considered.}

  \item{modelList}{a list contained 3 models fits using \code{\link{drm}} with the model fit
  for single mixture ratio being the first element, followed by the 2 model fits of the pure
  substances.}
  
  \item{EDvec}{ a vector of numeric values between 0 and 100 (percentages) coresponding to the
  effect levels of interest.}

  \item{EDonly}{a logical value indicating whether or not only combination indices for 
  effective doses should be calculated.}
  
  \item{effList}{a list returned by \code{\link{CIcompX}}.}
  
  \item{indAxis}{a character indicating whether effective doses ("ED") or effects
  ("EF") should be plotted.}
  
  \item{caRef}{a logical value indicating whether or not a reference line for concentration
  addition should be drawn.}
  
  \item{showPoints}{A logical value indicating whether or not estimated combination indices
  should be plotted.}
  
  \item{add}{a logical value specifying if the plot should be added to the existing plot.}
  
  \item{ylim}{a numeric vector of length 2 giving the range for the y axis.}
  
  \item{...}{additional graphical arguments.}
}

\details{
\code{\link{CIcomp}} calculates the classical combination index for effective doses whereas \code{\link{CIcompX}} calculates the combination index also for effects as proposed by
Martin-Betancor et al. (2015); for details and examples using "drc" see the supplementary material of this paper. The function \code{\link{plotFACI}} may be used to visualize the
calculated combination index as a function of the fraction affected.
}

\value{
\code{\link{CIcomp}} returns a matrix which one row per ED value. Columns contain
estimated combination indices, their standard errors and 95\% confidence intervals,
p-value for testing CI=1, estimated ED values for the mixture data and assuming 
concentration addition (CA) with corresponding standard errors.

\code{\link{CIcompX}} returns similar output both for effective doses and effects (as a
list of matrices). 
}

\references{
Martin-Betancor, K. and Ritz, C. and Fernandez-Pinas, F. and Leganes, F. and Rodea-Palomares, I. (2015) 
Defining an additivity framework for mixture research in inducible whole-cell biosensors,
\emph{Scientific Reports}
\bold{17200}.
}

\author{
Christian Ritz and Ismael Rodea-Palomares
}

%\note{}

\seealso{
See \code{\link{mixture}} for simultaneous modelling of several mixture ratios, but only at the ED50 level.

See also the help page for \code{\link{metals}}.
}

\examples{
## Fitting marginal models for the 2 pure substances
acidiq.0 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 0), fct = LL.4())
acidiq.100 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 100), fct = LL.4())

## Fitting model for single mixture with ratio 17:83
acidiq.17 <- drm(rgr ~ dose, data = subset(acidiq, pct == 17 | pct == 0), fct = LL.4())
    
## Calculation of combination indices based on ED10, ED20, ED50
CIcomp(0.17, list(acidiq.17, acidiq.0, acidiq.100), c(10, 20, 50))  
## CI>1 significantly for ED10 and ED20, but not so for ED50

}
\keyword{models}
\keyword{nonlinear}

\concept{antagonism mixture synergy}
