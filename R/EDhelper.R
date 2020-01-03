"EDhelper" <- function(parmVec, respl, reference, typeCalc, cond = TRUE)
{
    ## Works for log-logistic type dose-response models
  
    ## Converting absolute to relative
    if (typeCalc == "absolute") 
    {
        p <- 100 * ((parmVec[3] - respl) / (parmVec[3] - parmVec[2]))
#        typeCalc <- "relative"
    } else {  
        p <- respl
    }
    ## Swapping p for an increasing curve
    if (cond)
    {
        if ((typeCalc == "relative") && (parmVec[1] < 0) && (reference == "control"))
        {
            p <- 100 - p
        }
    } else {
        if ((typeCalc == "relative") && (reference == "control"))
        {
            p <- 100 - p
        }
    }
    p
}

## Used for FP models
"EDhelper2" <- function(parmVec, respl, reference, typeCalc, increasing)
{
  ## Converting absolute to relative
  if (typeCalc == "absolute") 
  {
    p <- 100 * (1 - (parmVec[3] - respl) / (parmVec[3] - parmVec[2]))

  } else {  

    if (increasing) {p <- respl} else {p <- 100 - respl}
  }

  p
}