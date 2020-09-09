#
# Wrap-up version of auc.roca
#
#
#

rm(list=ls())
source("FUN_GreyZone_v2.R")

auc.nocheck<- function(ss) {
    
    if (is.unsorted(ss[,1])) {
        ss[,1] <- rev(ss[,1])
        ss[,2] <- rev(ss[,2])
    }
    
    sp = ss[,1]
    se = ss[,2]
    diffs.x <- sp[-1] - sp[-length(sp)]
    means.vert <- (se[-1] + se[-length(se)])/2
    auc <- sum(means.vert * diffs.x)
    return(auc)
}

auc.wrapup <- function(controls, cases, direction = c("auto", "<",">")) {
    
    # set direction by median
    direction <- match.arg(direction)
    
    if (direction == "auto" && median(controls[,1]) <= median(cases[,1])) {
        direction <- "<"
    }
    else if (direction == "auto" && median(controls[,1]) > median(cases[,1])) {
        direction <- ">"
    } 
    
    
    # create legit threshold candidate
    tmat=thresholds.new(controls,cases)
    tt=expand.grid(tmat[,1],tmat[,2],tmat[,3])
    tt=tt[tt[,1]<=tt[,2],]
    
    ### single step
    st=tmat[,1]
    ss.single= t(sapply(st, perfs.single, controls = controls, cases = cases, direction = direction))
    
    ### GreyZone
    # sensitivity and specificity
    kk=apply(tt, 1, perfs.new, controls=controls, cases=cases, direction=direction)
    
    # different type of curves (sensitivity and specificity)
    ss.list = curveFinder(t(kk)[,c(1:2)])

    # calculate auc
    auc.single = auc.nocheck(ss.single)
    auc.upper = auc.nocheck(ss.list$ss.upper)
    auc.lower = auc.nocheck(ss.list$ss.lower)
    
    return(c(auc.single,
             auc.upper,
             auc.lower))
}

