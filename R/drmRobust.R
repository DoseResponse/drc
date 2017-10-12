"drmRobust" <- function(robust, fctCall, lenData, lenPar)
{ 

    ## Finding robust scale estimate for trimmed and winsorised means and tukey
    if (robust%in%c("trimmed", "tukey", "winsor"))
    {
        call1 <- fctCall
        call1$"robust" <- "lts"

        scaleEst <- mad(residuals(eval(call1, parent.frame())), 0)
    }
    

    ## Defining distance functions
    quadratic <- function(x) {x*x}
    
    lms <- function(x) {median(x*x)}

    noRes <- floor((lenData+lenPar+1)/2)
    lts <- function(x) {sum(((x[order(x)])[1:noRes])^2)}
 
    metricTrim <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        cVal <- 1.345
        retVal <- x*x
        
#        indexVec <- abs(x) > cVal
#        print(x)
#        sumVec <- sum(indexVec)
#        print(sumVec)
#        
#        if (sumVec>0) {retVal[indexVec] <- rep(cVal * cVal, sumVec)}
        retVal[abs(x) > cVal] <- cVal^2        
        retVal
    }

    metricWinsor <- function(x) 
    {       
        if (all(is.na(x))) {return(x)}
    
        cVal <- 1.345

#        print(mad(x,0))
#        scaleEst <- (median(abs(x))/0.6745)  # overrules general scale estimate

#        scaleEst <- 9.03  # 9.055 for phones data set

        x <- x/scaleEst

        retVal <- x*x

        indexVec <- abs(x) > cVal
#        sumVec <- sum(indexVec)
        
#        if (sumVec>0) {retVal[indexVec] <- (cVal * (2 * abs(x) - cVal))[indexVec]}
        retVal[indexVec] <- (cVal * (2 * abs(x) - cVal))[indexVec]  
#        retVal[abs(x)>c] <- (c*(2*abs(x)-c))[abs(x)>c] 
        retVal
    }

    tukeyBiweight <- function(x) 
    {
        if (all(is.na(x))) {return(x)}
    
        x <- x/scaleEst
    
        Rval <- 4.685
        retVal <- (x^6)/(Rval^4) - 3*(x^4)/(Rval^2) + 3*x*x
        
#        indexVec <- abs(x) > Rval
#        sumVec <- sum(indexVec)
#        
#        if (sumVec>0) {retVal[indexVec] <- rep(Rval * Rval, sumVec)}
        retVal[abs(x) > Rval] <- Rval * Rval
        retVal
    }


    ## Assigning objective function
    robustFct <- switch(robust, 
                        mean = quadratic, 
                        median = abs, 
                        trimmed = metricTrim, 
                        tukey = tukeyBiweight, 
                        winsor = metricWinsor, 
                        lms = lms, 
                        lts = lts)
    
    return(robustFct)    
}
