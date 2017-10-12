backfit <- function(drcObject)
{
    DL <- drcObject$dataList
    DLdose <- DL$dose
    meansVec <- tapply(DL$origResp, DLdose, mean, na.rm = TRUE) 
    # arranged according to ascending dose values
    # therefore unique doses are sorted below

    backfitValues <- ED(drcObject, meansVec, type = "absolute", 
                        display = FALSE, multcomp = FALSE)[, 1, drop = FALSE]

#     colnames(backfitValues) <- "backfit"
#     rownames(backfitValues) <- sort(unique(DLdose))
#     backfitValues
    
    retMat <- cbind(dose = sort(unique(DLdose)), backfit = backfitValues)
    rownames(retMat) <- NULL
    return(retMat)
}