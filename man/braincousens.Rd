\name{braincousens}

\alias{braincousens}

\title{The Brain-Cousens hormesis models}

\description{
  'braincousens' provides a very general way of specifying Brain-Cousens' 
  modified log- logistic model for describing hormesis, under various constraints on the parameters.
}

\usage{
  braincousens(fixed = c(NA, NA, NA, NA, NA), 
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"), ssfct = NULL,
  fctName, fctText)
}

\arguments{
  \item{fixed}{numeric vector. Specifies which parameters are fixed and at what value they are fixed. 
  NAs for parameter that are not fixed.}
  \item{names}{a vector of character strings giving the names of the parameters (should not contain ":"). 
  The default is reasonable (see under 'Usage'). The order of the parameters is: b, c, d, e, f (see under 'Details').} 
  \item{method}{character string indicating the self starter function to use.}  
  \item{ssfct}{a self starter function to be used.}  
  \item{fctName}{optional character string used internally by convenience functions.}  
  \item{fctText}{optional character string used internally by convenience functions.}                            
}

\details{
  The Brain-Cousens model is given by the expression
  \deqn{ f(x) = c + \frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}}
  which is a five-parameter model.
  
  It is a modification of the four-parameter logistic curve to take hormesis into account proposed
  by Brain and Cousens (1989).
}

\value{
  The value returned is a list containing the non-linear function, the self starter function, 
  the parameter names and additional model specific objects.  
}

\references{
  Brain, P. and Cousens, R. (1989) An equation to describe dose responses 
  where there is stimulation of growth at low doses, 
  \emph{Weed Research}, \bold{29}, 93--96.
}
            
\author{Christian Ritz}

\note{
  This function is for use with the function \code{\link{drm}}.

  The convenience functions of \code{braincousens} are \code{\link{BC.4}} and \code{\link{BC.5}}. These functions
  should be used rather than \code{braincousens} directly.
}

%\seealso{}

%\examples{}

\keyword{models}
\keyword{nonlinear}

\concept{hormesis hormetic effect initial stimulation}
