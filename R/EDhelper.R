"EDhelper" <- function(parmVec, respl, reference, typeCalc, cond = TRUE)
# cond = TRUE for log-logistic and Weibull type 1 models
# cond = FALSE for log-normal and Weibull type 2 models
{
    ## Converting absolute to relative
    if (typeCalc == "absolute") {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
        if (p < 0 || p > 100) {
            warning("The specified response value is outside the range of the model fit")
        }
    } else {  
        p <- respl
    }

    ## Swapping p for an increasing fitted dose-response curve
    if (cond) {
        if ((typeCalc == "relative") && (parmVec[1] < 0) && (reference == "control")) {
            p <- 100 - p
        }
    } else {
        if ((typeCalc == "relative") && (parmVec[1] > 0) && (reference == "control")) {
            p <- 100 - p
        }
    }
    p
}

## Used for FP models
"EDhelper2" <- function(parmVec, respl, reference, typeCalc, increasing)
{
  ## Converting absolute to relative
  if (typeCalc == "absolute") {
    p <- 100 * (1 - (parmVec[3] - respl) / (parmVec[3] - parmVec[2]))

  } else {
    if (increasing) {p <- respl} else {p <- 100 - respl}
  }
  p
}