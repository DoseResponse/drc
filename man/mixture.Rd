\name{mixture}

\alias{mixture}

\title{Fitting binary mixture models}

\description{
  'mixture' fits a concentration addition, Hewlett or Voelund model to data from binary mixture toxicity experiments.
}

\usage{
  mixture(object, model = c("CA", "Hewlett", "Voelund"), start, startm, control = drmc())
}

\arguments{
  \item{object}{object of class 'drc' corresponding to the model with freely varying EC50 values.}
  \item{model}{character string. It can be "CA", "Hewlett" or "Voelund".}
  \item{start}{optional numeric vector supplying starting values for all parameters in the mixture model.}
  \item{startm}{optional numeric vector supplying the lambda parameter in the Hewlett model or 
  the eta parameters (two parameters) in the Voelund model.}  
  \item{control}{list of arguments controlling constrained optimisation (zero as boundary), 
  maximum number of iteration in the optimisation,
  relative tolerance in the optimisation, warnings issued during the optimisation.}
}

\details{
  The function is a wrapper to \code{\link{drm}}, implementing the models described in Soerensen et al. (2007). 
  See the paper for a discussion of the merits of the different models.
  
  Currently only the log-logistic models are available. Application of Box-Cox transformation is not yet available.
}

\value{
  An object of class 'drc' with a few additional components. 
}

\references{
  Ritz, C. and Streibig, J. C. (2014)
  From additivity to synergism - A modelling perspective
  \emph{Synergy}, \bold{1}, 22--29.
}

\author{Christian Ritz}

\seealso{
  The examples in \code{\link{acidiq}} (the Hewlett model), \code{\link{glymet}} (dose/concentration addition) 
  and \code{\link{mecter}} (the Voelund model).
}

%\examples{}

\keyword{models}
\keyword{nonlinear}
