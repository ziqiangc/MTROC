kde.wrapup <- function(x, supp, xmin, xmax) {
    H <- Hpi.diag.kcde(x = x)
    Fhat <- kde(x, H=H*0.7, supp = supp, xmin = xmin, xmax = xmax)
    return(Fhat)
}

kcde.wrapup <- function(Fhat, eval.points) {
    diffe1 <- abs(diff(Fhat$eval.points[[1]]))
    diffe2 <- abs(diff(Fhat$eval.points[[2]]))
    
    Fhatsum <- matrix(apply(Fhat$estimate, 1, sum), ncol = ncol(Fhat$estimate), 
                      nrow = nrow(Fhat$estimate), byrow = TRUE)
    Fhat$estimate <- (Fhatsum - apply(Fhat$estimate, 
                                      1, cumsum)) * c(diffe1[1], diffe1)
    Fhatsum <- matrix(apply(Fhat$estimate, 1, sum), ncol = ncol(Fhat$estimate), 
                      nrow = nrow(Fhat$estimate), byrow = TRUE)
    Fhat$estimate <- (Fhatsum - apply(t(Fhat$estimate), 
                                      2, cumsum)) * c(diffe2[1], diffe2)

    Fhat$estimate <- Fhat$estimate/max(Fhat$estimate)
    
    Fhat$estimate <- predict(Fhat, x = eval.points)
    Fhat$eval.points <- eval.points
    
    return(Fhat)
}



