\name{leaflength}

\alias{leaflength}

\docType{data}

\title{Leaf length of barley}

\description{
  In an experiment barley was grown in a hydroponic solution with a herbicide.
}

\usage{data(leaflength)}

\format{
  A data frame with 42 observations on the following 2 variables.
  \describe{
    \item{\code{Dose}}{a numeric vector}
    \item{\code{DW}}{a numeric vector}
  }
}

\details{
  The dataset exhibits a large hormetical effect.
}

\source{
  Nina Cedergreen, Royal Veterinary and Agricultural University, Denmark.
}

%\references{}

\examples{

## Fitting a hormesis model
leaflength.crs4c1 <- drm(DW ~ Dose, data = leaflength, fct = CRS.4c())
plot(fitted(leaflength.crs4c1), residuals(leaflength.crs4c1))

leaflength.crs4c2 <- boxcox(drm(DW ~ Dose, data = leaflength, fct = CRS.4c()), 
method = "anova", plotit = FALSE)
summary(leaflength.crs4c2)

## Plottinf fitted curve and original data
plot(leaflength.crs4c2, broken = TRUE, conLevel = 0.001, type = "all", legend = FALSE, 
ylab = "Produced leaf length (cm)", xlab = "Metsulfuron-methyl (mg/l)",
main = "Hormesis: leaf length of barley")

}
\keyword{datasets}
