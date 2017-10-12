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
  Specifying only a single object gives a test for lack-of-fit, comparing the non-linear regression
  model to a more general one-way or two-way ANOVA model.
  
  If two objects are specified a test for reduction from the larger to the smaller model is given. (This only makes statistical 
  sense if the models are nested, that is: one model is a submodel of the other model.)
}

\value{
  An object of class 'anova'.
}

\references{
  Bates, D. M. and Watts, D. G. (1988) 
  \emph{Nonlinear Regression Analysis and Its Applications}, 
  New York: Wiley \& Sons (pp. 103--104)
} 

\author{Christian Ritz}

\seealso{
  For comparison of nested or non-nested model the function \code{\link{mselect}}can also be used.
  
  The function \code{\link{anova.lm}} for linear models.
}

\examples{

## Comparing a Gompertz three- and four-parameter models using an F test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
anova(ryegrass.m2, ryegrass.m1)  # reduction to 'W1.3' not possible (highly significant)   

anova(ryegrass.m2, ryegrass.m1, details = FALSE)  # without details

}
\keyword{models}
\keyword{nonlinear}

\concept{extra sum-of-squares F-test}
\concept{approximate F-test}