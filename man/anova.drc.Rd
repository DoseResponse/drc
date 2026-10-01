\name{anova.drc}

\alias{anova.drc}

\title{ANOVA for dose-response model fits}

\description{
  'anova' produces an analysis of variance table for one or two non-linear model fits.
}

\usage{

  \method{anova}{drc}(object, ..., details = TRUE, test = NULL)
}

\arguments{
  \item{object}{an object of class 'drc'.}
  \item{...}{additional arguments.}
  \item{details}{logical indicating whether or not details on the models compared should be displayed. 
  Default is TRUE (details are displayed).}
  \item{test}{a character string specifying the test statistic to be applied.
              Use "od" to assess overdispersion for binomial data.}
}

\details{
  If a single object is specified a message is suggesting to use the function \code{\link{modelFit}} for test of lack-of-fit, 
  comparing the dose-response model to a more general one-way or two-way ANOVA model.
  
  If two objects are specified a test for reduction from the larger (more general) model to the smaller, reduced 
  model is given. This only makes sense if the models are nested: the reduced model is a submodel of the larger model.
}

\value{
  An object of class 'anova'.
}

\references{
  Bates, D. M. and Watts, D. G. (1988) 
  \emph{Nonlinear Regression Analysis and Its Applications}, 
  New York: Wiley and Sons (pp. 103--104)
} 

\author{Christian Ritz}

\seealso{
  For comparison of nested or non-nested model the function \code{\link{mselect}} can also be used.
  
  The function \code{\link{anova.lm}} for linear models.
}

\examples{

## Comparing a four-parameter Weibull type 1 model 
##  to a reduced three-parameter Weibull type 1 using an F test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
anova(ryegrass.m2, ryegrass.m1)  
# reduction to 'W1.3' not possible (highly significant)   

anova(ryegrass.m2, ryegrass.m1, details = FALSE)  
# not showing details

}
\keyword{models}
\keyword{nonlinear}

\concept{extra sum-of-squares F-test}
\concept{approximate F-test}