"rse" <- function(object, resvar = FALSE)
{
    if (!is.null(object$"objList"))
    {
        fitValue <- object$"minval"
    } else {
        fitValue <- object$"fit"$"value"
    }

    rse <- switch(object$"type",
    "continuous" = fitValue / df.residual(object),
    "binomial" = NA,
    "Poisson" = NA,
    "negbin1" = NA,
    "negbin2" = NA,
    "event" = NA,
    "ssd" = NA,
    "standard" = fitValue / df.residual(object))
    
    if (resvar) 
    {
        rse
    } else {
        sqrt(rse)
    }
}