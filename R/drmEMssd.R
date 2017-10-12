"drmEMssd" <- 
function(dose, resp, multCurves, startVec, weightsVec, doseScaling = 1, multCurves2 = NULL)
{
    censYes <- (ncol(as.matrix(dose)) > 1)  #(!is.null(ncol(dose))) & (ncol(dose == 2))
    if (censYes)
    {
        dose1 <- dose[, 1]
        dose2 <- dose[, 2]      
        notCens <- dose1 == dose2
    }
  
    ## Defining the objective function  
    opfct <- function(cVal)
    {
#        -sum(log(multCurves(dose / doseScaling, cVal)))
        # not using resp   
 
      # Handling censoring 
      if (censYes)
      {
        fValues <- multCurves(dose1 / doseScaling, cVal)[notCens]
        Fvalues1 <- multCurves2(dose1 / doseScaling, cVal)[!notCens]
        Fvalues2 <- multCurves2(dose2 / doseScaling, cVal)[!notCens]
        #print(multCurves(dose1 / doseScaling, cVal))
#        print(fValues)
#        print(Fvalues1)
#        print(Fvalues2)
        -sum(log(fValues)) + (-sum(log(Fvalues2 - Fvalues1)))
      } else {
        -sum(log(multCurves(dose / doseScaling, cVal)))
        
      } 
    }
 

    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
        c(-object$"fit"$value, object$"sumList"$"df.residual")
    }
    
       
    ## Defining functions returning the residual variance, the variance-covariance and the fixed effects estimates
    rvfct <- NULL

    vcovfct <- function(object)
    {
        solve(object$fit$hessian)    
    }
    
    parmfct <- function(fit, fixed = TRUE)
    {
        fit$par
    }


    ## Returning list of functions
    return(list(llfct = llfct, opfct = opfct, ssfct = ssfct, rvfct = rvfct, vcovfct = vcovfct, 
    parmfct = parmfct))
}


"drmLOFssd" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
