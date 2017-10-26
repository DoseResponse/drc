"drmEMnegbin" <- 
function(dose, resp, multCurves, startVec, weightsVec, doseScaling = 1, dist.type = 1)
{

    ## Finding indices for doses that give contribution to likelihood function
#    iv <- ( (multCurves(dose, startVec) > zeroTol) & (multCurves(dose, startVec) < 1-zeroTol) )


    ## Defining the objective function  
    if (dist.type == 1)
    {  
        opfct <- function(cVal)
        {
            sizeVal <- tail(cVal, 1)
            pVal <- 1 / (1 + weightsVec * multCurves(dose / doseScaling, head(cVal, -1)) * exp(sizeVal)) 
#            print(c(sizeVal, pVal, -sum(dnbinom(resp, exp(-sizeVal), pVal, log = TRUE))))
            -sum(dnbinom(resp, exp(-sizeVal), pVal, log = TRUE))
        }
    }

    if (dist.type == 2)
    {  
        opfct <- function(cVal)
        {
            sizeVal <- tail(cVal, 1)
            pVal <- 1 / (1 + weightsVec * exp(sizeVal)) 
            -sum(dnbinom(resp, exp(-sizeVal) * multCurves(dose / doseScaling, head(cVal, -1)), 
                         pVal, log = TRUE))
        }
    }
  
    
    ## Defining self starter function
    ssfct <- NULL


    ## Defining the log likelihood function
    llfct <- function(object)
    {
#        total <- (object$"data")[iv, 5]
#        success <- total*(object$"data")[iv, 2]    
#        c( sum(log(choose(total, success))) - object$"fit"$"ofvalue", object$"sumList"$"df.residual" )
        
        c(-object$"fit"$value, object$"sumList"$"df.residual"
        )  # adding scale constant
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


"drmLOFnegbin" <- function()
{
    return(list(anovaTest = NULL, gofTest = NULL))
}
