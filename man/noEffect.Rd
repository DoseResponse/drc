\name{noEffect}

\alias{noEffect}

\title{
  Testing if there is a dose effect at all
}

\description{
  A significance test is provided for the comparison of the dose-response model considered and the simple linear regression
  model with slope 0 (a horizontal regression line corresponding to no dose effect)
}

\usage{
noEffect(object)
}

\arguments{
  \item{object}{an object of class 'drc'.}
}

\details{
  Perhaps useful for screening purposes.
}

\value{
  The likelihood ratio test statistic and the corresponding degrees of freedom and p-value are reported.
}

%\references{
%}

\author{Christian Ritz}

%\note{
%}

\examples{

ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

noEffect(ryegrass.LL.4)
# p-value < 0.0001: there is a highly significant dose effect!

}

\keyword{models}
\keyword{nonlinear}

