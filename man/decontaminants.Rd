\name{decontaminants}

\alias{decontaminants}

\docType{data}

\title{
Performance of decontaminants used in the culturing of a micro-organism
}

\description{
The two decontaminants 1-hexadecylpyridium chloride and oxalic acid were used. Additionally there was a control group (coded as concentration 0 and only included under oxalic acid).
}

\usage{data("decontaminants")}
\format{
  A data frame with 128 observations on the following 3 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of percentage weight per volume}
    \item{\code{count}}{a numeric vector of numbers of M. bovis colonies at stationarity}
    \item{\code{group}}{a factor with levels \code{hpc} and \code{oxalic} of the decontaminants used}
  }
}

\details{
These data examplify Wadley's problem: counts where the maximum number is not known. The data were analyzed by Trajstman (1989) using a three-parameter logistic model and then re-analyzed by Morgan and Smith (1992) using a three-parameter Weibull type II model. In both cases the authors adjusted for overdispersion (in different ways).

It seems that Morgan and Smith (1992) fitted separate models for the two decontaminants and using the control group for both model fits. In the example below a joint model is fitted where the control group is used once to determine a shared upper limit at concentration 0.
}

\source{
Trajstman, A. C. (1989) Indices for Comparing Decontaminants when Data Come from Dose-Response Survival and Contamination Experiments, \emph{Applied Statistics}, \bold{38}, 481--494. 
}

\references{
Morgan, B. J. T. and Smith, D. M. (1992) A Note on Wadley's Problem with Overdispersion, \emph{Applied Statistics}, \bold{41}, 349--354.
}

\examples{

## Wadley's problem using a three-parameter log-logistic model
decon.LL.3.1 <- drm(count~conc, group, data = decontaminants, fct = LL.3(), 
type = "Poisson", pmodels = list(~group, ~1, ~group))

summary(decon.LL.3.1)

plot(decon.LL.3.1)


## Same model fit in another parameterization (no intercepts)
decon.LL.3.2 <- drm(count~conc, group, data = decontaminants, fct=LL.3(), 
type = "Poisson", pmodels = list(~group-1, ~1, ~group-1))

summary(decon.LL.3.2)

}
\keyword{datasets}
