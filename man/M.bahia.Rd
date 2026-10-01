\name{M.bahia}

\alias{M.bahia}

\docType{data}

\title{Effect of an effluent on the growth of mysid shrimp}

\description{
  Juvenile mysid shrimp (\emph{Mysidopsis bahia}) were exposed to up to 32\% effluent
  in a 7-day survival and growth test. The average weight per treatment replicate of
  surviving organisms was measured.
}

\usage{data(M.bahia)}

\format{
  A data frame with 40 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of effluent concentrations (\%)}
    \item{\code{dryweight}}{a numeric vector of average dry weights (mg)}
  }
}

\details{
  The data are analysed in Bruce and Versteeg (1992) using a log-normal 
  dose-response model (using the logarithm with base 10).
  
  At 32\% there was complete mortality, and this justifies using a model where a lower asymptote
  of 0 is assumed.
}

\source{
  Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for modeling continuous toxicity data,
  \emph{Environ. Toxicol. Chem.}, \bold{11}, 1485--1494.
}

%\references{}

\examples{

M.bahia.m1 <- drm(dryweight~conc, data=M.bahia, fct=LN.3())

## Variation increasing
plot(fitted(M.bahia.m1), residuals(M.bahia.m1))

## Using transform-both-sides approach
M.bahia.m2 <- boxcox(M.bahia.m1, method = "anova")
summary(M.bahia.m2)  # logarithm transformation

## Variation roughly constant, but still not a great fit
plot(fitted(M.bahia.m2), residuals(M.bahia.m2))

## Visual comparison of fits
plot(M.bahia.m1, type="all", broken=TRUE)
plot(M.bahia.m2, add=TRUE, type="none", broken=TRUE, lty=2)

ED(M.bahia.m2, c(10,20,50), ci="fls")

## A better fit
M.bahia.m3 <- boxcox(update(M.bahia.m1, fct = LN.4()), method = "anova")
#plot(fitted(M.bahia.m3), residuals(M.bahia.m3))
plot(M.bahia.m3, add=TRUE, type="none", broken=TRUE, lty=3, col=2)
ED(M.bahia.m3, c(10,20,50), ci="fls")

}
\keyword{datasets}
