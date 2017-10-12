\name{getInitial}

\Rdversion{1.1}

\alias{getInitial}

\title{
  Showing starting values used
}

\description{
  Function for showing the starting values of the model parameters used when fitting a dose-response model
}

\usage{
  getInitial(object)
}

\arguments{
  \item{object}{object of class 'drc'}
}

%\details{}

\value{
  A vector of starting values for the model parameters used to initialize the estimation procedure.
}

%\references{}

\author{Christian Ritz}

\note{This function is masking the standard function in the stats package.}

%\seealso{}

%\examples{}

\keyword{models}
\keyword{nonlinear}