\name{vinclozolin}

\alias{vinclozolin}

\docType{data}

\title{Vinclozolin from AR in vitro assay}

\description{
  Dose-response experiment with vinclozolin in an AR reporter gene assay
}

\usage{data(vinclozolin)}

\format{
  A data frame with 53 observations on the following 3 variables.
  \describe{
    \item{\code{exper}}{a factor with levels \code{10509} \code{10821} \code{10828} \code{10904} \code{11023} \code{11106}}
    \item{\code{conc}}{a numeric vector of concentrations of vinclozolin}
    \item{\code{effect}}{a numeric vector of luminescense effects}
  }
}

\details{
  The basic dose-response experiment was repeated 6 times on different days. Chinese Hamster Ovary cells 
  were exposed to various concentrations of vinclozolin for 22 hours and the resulting luminescense effects
  were recorded.
  
  Data are part of mixture experiment reported in Nellemann \emph{et al} (2003).  
}

\source{
  Nellemann C., Dalgaard M., Lam H.R. and Vinggaard A.M. (2003)
  The combined effects of vinclozolin and procymidone do not deviate from expected additivity \emph{in vitro} 
  and \emph{in vivo}, \emph{Toxicological Sciences}, \bold{71}, 251--262.
}

%\references{}

\examples{


vinclozolin.m1 <- drm(effect~conc, exper, data=vinclozolin, fct = LL.3())
plot(vinclozolin.m1, xlim=c(0,50), ylim=c(0,2800), conLevel=1e-4)

vinclozolin.m2 <- drm(effect~conc, data=vinclozolin, fct = LL.3())
plot(vinclozolin.m2, xlim=c(0,50), conLevel=1e-4, add=TRUE, type="none", col="red")

## Are the ED50 values indetical across experiments?
vinclozolin.m3 <- update(vinclozolin.m1, pmodels=data.frame(exper, exper, 1))
anova(vinclozolin.m3, vinclozolin.m1)  # No!

}
\keyword{datasets}
