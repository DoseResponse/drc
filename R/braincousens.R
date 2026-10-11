"braincousens" <- function(
    fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), 
    method = c("1", "2", "3", "4"), ssfct = NULL,
    fctName, fctText) {

    ## Checking arguments
    numParm <- 5
    if (!is.character(names) | !(length(names)==numParm)) {stop("Not correct 'names' argument")}
    if (!(length(fixed)==numParm)) {stop("Not correct 'fixed' argument")}    

#    if (!is.logical(useDer)) {stop("Not logical useDer argument")}
#    if (useDer) {stop("Derivatives not available")}

    notFixed <- is.na(fixed)
    parmVec <- rep(0, numParm)
    parmVec[!notFixed] <- fixed[!notFixed]
    parmVec1 <- parmVec
    parmVec2 <- parmVec
    
    ## Defining the non-linear function
    fct <- function(dose, parm) 
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm
        
        parmMat[,2]+(parmMat[,3]+parmMat[,5]*dose-parmMat[,2])/(1+exp(parmMat[,1]*(log(dose)-log(parmMat[,4]))))
    }


#    ## Defining value for control measurements (dose=0)
#    confct <- function(drcSign)
#    {
#        if (drcSign>0) {conPos <- 2} else {conPos <- 3}
#        confct2 <- function(parm)
#        { 
#            parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
#            parmMat[, notFixed] <- parm
#            parmMat[, conPos]
#        }
#        return(list(pos=conPos, fct=confct2))
#    } 
#
#
#    ## Defining flag to indicate if more general ANOVA model
##    anovaYes <- TRUE
#    binVar <- all(fixed[c(2, 3, 5)]==c(0, 1, 1))
#    if (is.na(binVar)) {binVar <- FALSE}
#    if (!binVar) {binVar <- NULL}
#    anovaYes <- list(bin = binVar, cont = TRUE)


    ## Defining the self starter function
    if (!is.null(ssfct))
    {
        ssfct <- ssfct
    } else {
#        ssfct <- braincousens.ssf(method, fixed)
        ssfct <- function(dframe)
        {
            initval <- llogistic(method = method)$ssfct(dframe)   
            initval[5] <- 0
    
            return(initval[notFixed])
        }           
    }


#    ## Defining parameter to be scaled
#    if ( (scaleDose) && (is.na(fixed[4])) ) 
#    {
#        scaleInd <- sum(is.na(fixed[1:4]))
#    } else {
#        scaleInd <- NULL
#    }
   
   
    ## Defining names
    names <- names[notFixed]


    ## Defining derivatives
    deriv1 <- function(dose, parm)
    {
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow=TRUE)
        parmMat[, notFixed] <- parm

        t1 <- parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*dose
        t2 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t3 <- 1 + t2                          
        t4 <- (1 + t2)^(-2)

        cbind( -t1*xlogx(dose/parmMat[, 4], parmMat[, 1])*t4, 
               1 - 1/t3, 
               1/t3, 
               t1*t2*(parmMat[, 1]/parmMat[, 4])*t4, 
               dose/t3 )[, notFixed]
    }
        
    deriv2 <- NULL


    ## Defining the first derivative (in the dose)
    derivx <- function(dose, parm) {
  
        parmMat <- matrix(parmVec, nrow(parm), numParm, byrow = TRUE)
        parmMat[, notFixed] <- parm
  
        t1 <- exp(parmMat[, 1]*(log(dose) - log(parmMat[, 4])))
        t2 <- 1 + t1
  
        parmMat[, 5] / t2 - (parmMat[, 3] - parmMat[, 2] + parmMat[, 5]*dose) * t1 * parmMat[, 1] / (t2^2 * dose)
    }


    ## Defining parameter check function
    checkParm <- function(parmV){
        if (parmV[1] <= 1) {stop("Brain-Cousens model with b<=1 not meaningful")}
        if (parmV[5] < 0) {stop("Brain-Cousens model with f<0 not meaningful")}
    }


    ## Defining the ED function
    edfct <- function(parm, respl, reference, type, lower = NULL, upper = NULL, ...) {
    #  lower = 1e-3, upper = 1000, ...)
    # lower is used for finding the maximum hormesis, upper is used for finding the EDp
        parmVec[notFixed] <- parm
        checkParm(parmVec)

        p <- EDhelper(parmVec, respl, reference, type)
        tempVal <- 1 - p/100  # (100-p)/100

        ## Defining the numerator of derivative of the BC model function wrt. dose
        derx <- function(dose) {
            expVal <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            parmVec[5]*(1+expVal*(1-parmVec[1]))-(parmVec[3]-parmVec[2]+parmVec[5]*dose)*expVal*parmVec[1]/dose
        }
        zepVal <- ((parmVec[5]/(parmVec[3]-parmVec[2]))*parmVec[4]^(parmVec[1]))^(1/(parmVec[1]-1))
         
        if (is.null(lower)) {
            xvalues <- 10^(runif(100, -10, log10(zepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec<0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        }
        #interval <- c(lower, upper)     
        maxAt <- uniroot(derx, lower = lower, upper = zepVal)$root   
        if (is.null(upper)) {
            upper <- parmVec[4] * 10000
        }
        eqn <- function(dose) {tempVal*(1+exp(parmVec[1]*(log(dose)-log(parmVec[4]))))-(1+parmVec[5]*dose/(parmVec[3]-parmVec[2]))}
        
        EDp <- uniroot(eqn, lower = maxAt, upper = upper)$root

        tempVal1 <- exp(parmVec[1] * (log(EDp) - log(parmVec[4])))
        tempVal2 <- parmVec[3] - parmVec[2]
        derNum <- c(-tempVal * tempVal1 * (log(EDp) - log(parmVec[4])), 
                     ((1 - tempVal) * tempVal2 + parmVec[5] * EDp) / (tempVal * tempVal2^2),
                     1 / tempVal2, 
                     tempVal * tempVal1 * parmVec[1] / parmVec[4],
                     EDp / tempVal2)
        derDenom <- tempVal * tempVal1 *parmVec[1] / EDp - parmVec[5] / tempVal2
        EDder <- derNum / derDenom
        
        return(list(EDp, EDder[notFixed]))
    }


    ## Finding the maximal hormesis
    maxfct <- function(parm, lower = NULL, noRan = 100) {
        # Note that the code is mostly copied from the edfct function above, 
        #   but with some changes to find the maximum hormesis

        parmVec[notFixed] <- parm
        checkParm(parmVec)

        ## Defining the numerator of derivative of the BC model function wrt. dose
        derx <- function(dose) {
            expVal <- exp(parmVec[1]*(log(dose)-log(parmVec[4])))
            parmVec[5]*(1+expVal*(1-parmVec[1]))-(parmVec[3]-parmVec[2]+parmVec[5]*dose)*expVal*parmVec[1]/dose
        }
        
        zepVal <- ((parmVec[5]/(parmVec[3]-parmVec[2]))*parmVec[4]^(parmVec[1]))^(1/(parmVec[1]-1))
        if (is.null(lower)) {
            xvalues <- 10^(runif(noRan, -10, log10(zepVal)))
            dffVec <- derx(xvalues)
            dffVec[dffVec<0] <- NA
            lower <- xvalues[which.min(abs(dffVec))]
        }  
        maxAt <- uniroot(derx, lower = lower, upper = zepVal)$root

        return(c(maxAt, fct(maxAt, matrix(parm, 1, length(names)))))
    }


    returnList <- list(
        fct = fct, ssfct = ssfct, names = names, 
        deriv1 = deriv1, deriv2 = deriv2, derivx = derivx, 
        edfct = edfct, maxfct = maxfct,
        name = ifelse(missing(fctName), as.character(match.call()[[1]]), fctName),
        text = ifelse(missing(fctText), "Brain-Cousens (hormesis)", fctText),
        noParm = sum(is.na(fixed)), fixed = fixed)

    class(returnList) <- "braincousens"
    invisible(returnList)
}


"BC.4" <- function(
fixed = c(NA, NA, NA, NA), names = c("b", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 4)) {stop("Not correct 'names' argument")}
 
    return(braincousens(names=c(names[1], "c", names[2:4]), fixed = c(fixed[1], 0, fixed[2:4]),
    fctName = as.character(match.call()[[1]]),
    fctText = "Brain-Cousens (hormesis) with lower limit fixed at 0", ...))
}

"BC.5" <- function(
fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
{
    ## Checking arguments
    if (!is.character(names) | !(length(names) == 5)) {stop("Not correct 'names' argument")}

    return(braincousens(names = names, fixed = fixed,
    fctName = as.character(match.call()[[1]]), ...))
}

