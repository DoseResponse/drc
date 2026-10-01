\name{lettuce}

\alias{lettuce}

\docType{data}

\title{Hormesis in lettuce plants}

\description{
  Data are from an experiment where isobutylalcohol was dissolved in a nutrient solution in which lettuce
  (\emph{Lactuca sativa}) plants were grown. The plant biomass of the shoot was determined af 21 days.
}

\usage{data(lettuce)}

\format{
  A data frame with 14 observations on the following 2 variables.
  \describe{
    \item{conc}{a numeric vector of concentrations of isobutylalcohol (mg/l)}
    \item{weight}{a numeric vector of biomass of shoot (g)}
  }
}

\details{
  The data set illustrates hormesis, presence of a subtoxic stimulus at low concentrations.
}

\source{
  van Ewijk, P. H. and Hoekstra, J. A. (1993) 
  Calculation of the EC50 and its Confidence Interval When Subtoxic Stimulus Is Present,
  \emph{ECOTOXICOLOGY AND ENVIRONMENTAL SAFETY}, \bold{25}, 25--32.
}

\references{
  van Ewijk, P. H. and Hoekstra, J. A. (1994) 
  Curvature Measures and Confidence Intervals for the Linear Logistic Model,
  \emph{Appl. Statist.}, \bold{43}, 477--487.
}

\examples{

## Look at data
lettuce

## Monotonous dose-response model
lettuce.m1 <- drm(weight~conc, data=lettuce, fct=LL.3())

plot(lettuce.m1, broken = TRUE)

## Model fit in van Ewijk and Hoekstra (1994)
lettuce.m2 <- drm(weight~conc, data=lettuce, fct=BC.4())
modelFit(lettuce.m2)

plot(lettuce.m2, add = TRUE, broken = TRUE, type = "none", lty = 2)

## Hormesis effect only slightly significant
summary(lettuce.m2)

## Hormesis effect highly significant
##  compare with t-test for the "f" parameter in the summary output)
anova(lettuce.m1, lettuce.m2)

}
\keyword{datasets}
