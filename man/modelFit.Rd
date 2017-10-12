\name{modelFit}

\alias{modelFit}

\title{Assessing the model fit}

\description{
  Checking the fit of dose-response model by means of formal significance tests or graphical procedures.
}

\usage{
  modelFit(object, test = NULL, method = c("gof", "cum"))
}

\arguments{
  \item{object}{object of class 'drc'}
  \item{test}{character string defining the test method to apply}
  \item{method}{character string specifying the method to be used for assessing the model fit}
}

\details{
  Currently two methods are available. For continuous data the clasical lack-of-fit test is applied 
  (Bates and Watts, 1988). The test compares the dose-response model to a more general ANOVA model 
  using an approximate F-test. For quantal data the crude goodness-of-fit test based on Pearson's statistic is used. 
  
  None of these tests are very powerful. A significant test result is more alarming than a non-significant one. 
}

\value{
  An object of class 'anova' which will be displayed in much the same way as an ordinary ANOVA table.
}

\references{
  Bates, D. M. and Watts, D. G. (1988) 
  \emph{Nonlinear Regression Analysis and Its Applications}, 
  New York: Wiley \& Sons (pp. 103--104).
}

\author{Christian Ritz}

%\note{}

%\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }

\examples{

## Comparing the four-parameter log-logistic model 
##  to a one-way ANOVA model using an approximate F test
## in other words applying a lack-of-fit test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
modelFit(ryegrass.m1)   

}

\keyword{models}
\keyword{nonlinear}
