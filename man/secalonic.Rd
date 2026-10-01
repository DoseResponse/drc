\name{secalonic}

\alias{secalonic}

\docType{data}

\title{Root length measurements}

\description{
  Data stem from an experiment assessing the inhibitory effect of secalonic acids on plant growth.
}

\usage{data(secalonic)}

\format{
  A data frame with 7 observations on the following 2 variables:
  \describe{
    \item{\code{dose}}{a numeric vector containing dose values (mM)}
    \item{\code{rootl}}{a numeric vector containing root lengths (cm)}
  }
}

\details{
  For each dose the root length is an average three measurements.
}

\source{
  Gong, X. and Zeng, R. and Luo, S. and Yong, C. and Zheng, Q. (2004) Two new 
  secalonic acids from \emph{Aspergillus Japonicus} and their allelopathic effects on higher plants,
  \emph{Proceedings of International Symposium on Allelopathy Research and Application, 27-29 April, 
  Shanshui, Guangdong, China (Editors: R. Zeng and S. Luo)}, 209--217.
  
  Ritz, C (2009)
  Towards a unified approach to dose-response modeling in ecotoxicology
  \emph{To appear in Environ Toxicol Chem}.
}

%\references{}

\examples{

## Fitting a four-parameter log-logistic model
secalonic.m1 <- drm(rootl ~ dose, data = secalonic, fct = LL.4())
summary(secalonic.m1)

## Fitting a three-parameter log-logistic model
##  lower limit fixed at 0
secalonic.m2 <- drm(rootl ~ dose, data = secalonic, fct = LL.3())
summary(secalonic.m1)

## Comparing logistic and log-logistic models
## (Figure 1 in Ritz (2009))
secalonic.LL4 <- drm(rootl ~ dose, data = secalonic, fct = LL.4())
secalonic.L4 <- drm(rootl ~ dose, data = secalonic, fct = L.4())

plot(secalonic.LL4, broken=TRUE, ylim=c(0,7), xlab="Dose (mM)", ylab="Root length (cm)", 
cex=1.2, cex.axis=1.2, cex.lab=1.2, lwd=2)

plot(secalonic.L4, broken=TRUE, ylim=c(0,7), add=TRUE, type="none", lty=2, lwd=2)

abline(h=coef(secalonic.L4)[3], lty=3, lwd=2)

}
\keyword{datasets}
