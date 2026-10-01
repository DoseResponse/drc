\name{O.mykiss}

\alias{O.mykiss}

\docType{data}

\title{Test data from a 21 day fish test}

\description{
  Test data from a 21 day fish test following the guidelines OECD GL204, 
  using the test organism Rainbow trout \emph{Oncorhynchus mykiss}.
}

\usage{data(O.mykiss)}

\format{
  A data frame with 70 observations on the following 2 variables.
  \describe{
    \item{\code{conc}}{a numeric vector of concentrations (mg/l)}
    \item{\code{weight}}{a numeric vector of wet weights (g)}
  }
}

\details{
  Weights are measured after 28 days.
}

\source{
  Organisation for Economic Co-operation and Development (OECD) (2006)
  \emph{CURRENT APPROACHES IN THE STATISTICAL ANALYSIS OF ECOTOXICITY DATA: A GUIDANCE TO APPLICATION - ANNEXES},
  Paris (p. 65).
}

\references{
  Organisation for Economic Co-operation and Development (OECD) (2006)
  \emph{CURRENT APPROACHES IN THE STATISTICAL ANALYSIS OF ECOTOXICITY DATA: A GUIDANCE TO APPLICATION - ANNEXES},
  Paris (pp. 80--85).
}

\examples{

head(O.mykiss)

## Fitting exponential model
O.mykiss.m1 <- drm(weight ~ conc, data = O.mykiss, fct = EXD.2(), na.action = na.omit)
modelFit(O.mykiss.m1)
summary(O.mykiss.m1)

## Fitting same model with transform-both-sides approach
O.mykiss.m2 <- boxcox(O.mykiss.m1 , method = "anova")
summary(O.mykiss.m2)
# no need for a transformation

## Plotting the fit
plot(O.mykiss.m1, type = "all", xlim = c(0, 500), ylim = c(0,4),
xlab = "Concentration (mg/l)", ylab = "Weight (g)", broken = TRUE)

}
\keyword{datasets}
