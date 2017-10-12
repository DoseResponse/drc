## Calculating combination indices for x and y axes

CIcompX <- function(mixProp, modelList, EDvec, EDonly = FALSE)
{
    ## Checking the input
    if ( (mixProp < 0) | (mixProp > 1) ) {stop("Mixture proportion should be between 0 and 1")}
    if ( (!is.list(modelList)) | (length(modelList) != 3) ) {stop("Exactly 3 model fits should be provided in a list")}
    if (length(EDvec) < 1) {stop("At least effective dose level should be specified")} 

    ## Estimating effective doses and effects
    ese12Vec <- ED(modelList[[1]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]  # the mixture
    ese1Vec <- ED(modelList[[2]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]
    ese2Vec <- ED(modelList[[3]], EDvec, display = FALSE, bound = FALSE)[, 1:2, drop = FALSE]    
    eseMat <- as.matrix(cbind(ese12Vec[, 1], ese1Vec[, 1], ese2Vec[, 1], ese12Vec[, 2], ese1Vec[, 2], ese2Vec[, 2]))
    rownames(eseMat) <- as.character(EDvec)
    colnames(eseMat) <- c("ED.mix", "ED1", "ED2", "SE.mix", "SE1", "SE2")

    pred12 <- predict(modelList[[1]], data.frame(ese12Vec[, 1]), se.fit = TRUE)  
    # Note: Ignoring uncertainty in the ED values!
    pred1 <- predict(modelList[[2]], data.frame(ese1Vec[, 1]), se.fit = TRUE)
    pred2 <- predict(modelList[[3]], data.frame(ese2Vec[, 1]), se.fit = TRUE) 
    
    ## In case only a single ED level is specified 
    ## (as predict() then returns a vector, not a matrix)
    if (!is.matrix(pred12)) {pred12 <- matrix(pred12, nrow = 1)}
    if (!is.matrix(pred1)) {pred1 <- matrix(pred1, nrow = 1)}
    if (!is.matrix(pred2)) {pred2 <- matrix(pred2, nrow = 1)}
       
    predMat <- as.matrix(cbind(pred12[, 1], pred1[, 1], pred2[, 1], 
                               pred12[, 2], pred1[, 2], pred2[, 2]))
    rownames(predMat) <- as.character(EDvec)
    colnames(predMat) <- c("E.mix", "E1", "E2", "SE.mix", "SE1", "SE2")

    ## Calculating combination index for effective doses
    xCAfct <- function(eseVec)
    {
        combInd <- (mixProp*eseVec[1]/eseVec[2]) + ((1-mixProp)*eseVec[1]/eseVec[3]) 
        caDiff <- combInd - 1

        derivFct <- function(ecVec)
        {
            derivComp <- c(
            mixProp / ecVec[2] + (1-mixProp) / ecVec[3], 
            -mixProp * ecVec[1] / (ecVec[2]^2), 
            -(1 - mixProp) * ecVec[1] / ((ecVec[3]^2)))
            
            derivComp
        }
        derivVec <- derivFct(eseVec[1:3])
        diagVec <- diag(eseVec[4:6]^2)
 #       seCI <- sqrt(as.vector((-derivVec[2:3]) %*% (diagVec[2:3, 2:3]) %*% (-derivVec[2:3])))
        seDiff <- sqrt(as.vector(derivVec %*% diagVec %*% derivVec))

        derivFct2 <- function(ecVec)
        {
            addEff <- (mixProp / ecVec[2] + (1-mixProp) / ecVec[3])^{-1}
            derivComp2 <- c(0, (addEff^2) * mixProp / (ecVec[2]^2), (addEff^2) * (1 - mixProp) / (ecVec[3]^2))
        
            c(addEff, derivComp2)
        }
        dfRes2 <- derivFct2(eseVec[1:3])
        derivVec2 <- dfRes2[2:4]
        seDiff2 <- sqrt(as.vector(derivVec2 %*% diagVec %*% derivVec2))
#        print(seDiff2)

        retVec <- c(combInd, seDiff, c(combInd - 1.96 * seDiff, combInd + 1.96 * seDiff), 
                    caDiff, 2 * (1 - pnorm(abs(caDiff / seDiff))), dfRes2[1], seDiff2)
        names(retVec) <- c("combInd", "SE",  
                           "lowCI", "highCI", 
                           "CAdiff", "CAdiffp", "PredAdd", "sePredAdd")
        retVec
    }

# CI = {1/ [(pi/EX1) + (1-pi)/Ex2]}/eseVec[1]
# CI = {1/ [(pi/eseVec[2]) + (1-pi)/eseVec[3]]}/eseVec[1]
# CI = {1/ [(pi*eseVec[1]/eseVec[2]) + (1-pi)*eseVec[1]/eseVec[3]]}    
    
    ## Calculating combination index for effects
    yCAfct <- function(eseVec)
    {
        denomi <- (mixProp*eseVec[1]/eseVec[2]) + ((1-mixProp)*eseVec[1]/eseVec[3])
        combInd <- 1 / denomi
        caDiff <- combInd - 1

        derivFct <- function(ecVec)
        {
            derivComp <- (-1/(denomi^2)) * c(
            mixProp / ecVec[2] + (1-mixProp) / ecVec[3] , 
            -mixProp * ecVec[1] / (ecVec[2]^2), 
            -(1 - mixProp) * ecVec[1] / (ecVec[3]^2))
        }
        derivVec <- derivFct(eseVec[1:3])
        diagVec <- diag(eseVec[4:6]^2)
 #       seCI <- sqrt(as.vector((-derivVec[2:3]) %*% (diagVec[2:3, 2:3]) %*% (-derivVec[2:3])))
        seDiff <- sqrt(as.vector(derivVec %*% diagVec %*% derivVec))

        derivFct2 <- function(ecVec)
        {
            addEff <- (mixProp / ecVec[2] + (1-mixProp) / ecVec[3])^{-1}
            derivComp2 <- c(0, (addEff^2) * mixProp / (ecVec[2]^2), (addEff^2) * (1 - mixProp) / (ecVec[3]^2))
        
            c(addEff, derivComp2)
        }
        dfRes2 <- derivFct2(eseVec[1:3])
        derivVec2 <- dfRes2[2:4]
        seDiff2 <- sqrt(as.vector(derivVec2 %*% diagVec %*% derivVec2))                
        
        retVec <- c(combInd, seDiff, 
                    c(combInd - 1.96 * seDiff, combInd + 1.96 * seDiff),
                    caDiff, 2 * (1 - pnorm(abs(caDiff / seDiff))), dfRes2[1], seDiff2)
        names(retVec) <- c("combInd", "SE",
                           "lowCI", "highCI", 
                           "CAdiff", "CAdiffp", "PredAdd", "sePredAdd")
        retVec
    }
    CAxMat <- t(apply(eseMat, 1, xCAfct))
    rownames(CAxMat) <- EDvec    
    CAyMat <- t(apply(predMat, 1, yCAfct))  # not yCAfct
    rownames(CAyMat) <- EDvec   
    
    if (EDonly)
    {
        list(Effx = eseMat, CAx = CAxMat, EDvec = EDvec)
    } else {
        list(Effx = eseMat, Effy = predMat, CAx = CAxMat, CAy = CAyMat, EDvec = EDvec)
    }
}

CIcomp <- function(mixProp, modelList, EDvec)
{
    resLst <- CIcompX(mixProp, modelList, EDvec, EDonly = FALSE)
    resMt <- matrix(NA, 5, 5)

    resMt[c(1,2,4), 1:2] <- matrix(resLst[["Effx"]][c(2, 3, 1, 5, 6, 4)], 3, 2)
    resMt[3, 1:2] <- resLst[["CAx"]][7:8]
    resMt[5, ] <- resLst[["CAx"]][c(1:2, 4:6)]
    
    resMt <- cbind(resLst[["CAx"]][, -5], resLst[["Effx"]][, c(1, 4)])
    colnames(resMt)[6:7] <- c("ED.CA", "SE.CA")

    # colnames(resMt) <- c("Est", "SE", "CIlow", "CIupp", "p-val")
    # rownames(resMt) <- c("ED.A", "ED.B", "ED.CApred", "ED.mix", "CombInd")    
    
    resMt
}


plotFACI <- function(effList, indAxis = c("ED", "EF"), caRef = TRUE, 
                     showPoints = FALSE, add = FALSE, ylim, ...)
{
    indAxis <- match.arg(indAxis)
#    indMat <- CIcompX(mixProp, modelList, faValues)
  
    faValues <- effList[["EDvec"]]
    minfa <- min(faValues)
    faValues[faValues < 0] <- -(100 - abs(faValues[faValues < 0])) 
    
#    if (indAxis == "x") {plotMat <- indMat[[1]]} else {plotMat <- indMat[[2]]}
    plotMat <- switch(indAxis, ED = effList[["CAx"]], EF = effList[["CAy"]])
    xVec <- as.numeric(rownames(plotMat))
    xVec[faValues < 0] <- rev(xVec[faValues < 0])
    xVec[faValues > 0] <- rev(xVec[faValues > 0])
    yVec <- plotMat[, 1]
    seVec <- plotMat[, 2]
    if (caRef) 
    {
        yLimits <- c(0, max(yVec))
    } else {
        yLimits <- range(yVec)
    }
    if (!missing(ylim))
    {
        yLimits <- ylim
    }
    
    if (!add)
    {
        plot(xVec, yVec, type = "l", xlim = c(minfa - 1, max(faValues) + 1), ylim = yLimits, 
        xlab = "Fraction affected", ylab = "Combination index", ...)
        
        abline(h = 1, lty = 2, lwd = 2)
    
    } else {
        lines(xVec, yVec, ...)   
    }
    dispersion(xVec, yVec, ulim = 1.96 * seVec, arrow.gap = 0.15, ...)
    if (showPoints) 
    {
        points(xVec, yVec, pch = 1)
    }
    
    invisible(plotMat)
}
