\name{hatvalues.drc}

\alias{hatvalues.drc}
\alias{cooks.distance.drc}

\title{
  Model diagnostics for nonlinear dose-response models
}

\description{
   Hat values (leverage values) and Cook's distance are provided for nonlinear dose-response model fits using the same formulas
   as in linear regression but based on the corresponding but approximate quantities available for nonlinear models.   
}

\usage{
  \method{cooks.distance}{drc}(model, ...)

  \method{hatvalues}{drc}(model, ...)
}

\arguments{
  \item{model}{an object of class 'drc'.}
  \item{\dots}{additional arguments (not used).}
}

\details{
  Hat values and Cook's distance are calculated using the formula given by Cook et al. (1986) and McCullagh and Nelder (1989). 
  
  The output values can be assessed in the same way as in linear regression.
}

\value{
  A vector of leverage values (hat values) or values of Cook's distance (one value per observation).
}
\references{

  Cook, R. D. and Tsai, C.-L. and Wei, B. C. (1986)
  Bias in Nonlinear Regression,
  \emph{Biometrika}
  \bold{73}, 615--623.
  
  McCullagh, P. and Nelder, J. A. (1989)
  emph{Generalized Linear Models}, 
  Second edition, Chapman \& Hall/CRC.  
}

\author{
  Christian Ritz
}

%\note{
%}

\examples{

ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

hatvalues(ryegrass.LL.4)

cooks.distance(ryegrass.LL.4)

}
\keyword{models}
\keyword{nonlinear}