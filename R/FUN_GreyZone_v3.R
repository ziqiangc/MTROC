#
#
# new empirical function
#
#

thresholds.new=function (controls, cases)
{
    predictor1 <- c(controls[,1], cases[,1])
    predictor2 <- c(controls[,2], cases[,2])
    thresholds1 <- sort(unique(predictor1))
    thresholds1 <- (c(-Inf, thresholds1) + c(thresholds1, +Inf))/2
    theta1.lower <- thresholds1
    theta1.upper <- thresholds1
    thresholds2 <- sort(unique(predictor2))
    thresholds2 <- (c(-Inf, thresholds2) + c(thresholds2, +Inf))/2
    theta2 <- thresholds2
    return(matrix(c(theta1.lower,theta1.upper,theta2), ncol=3))
}


## single step
perfs.single=function (threshold, controls, cases, direction) 
{
    controls = controls[,1]
    cases = cases[,1]
    
    if (direction == ">") {
        tp <- sum(cases <= threshold)
        tn <- sum(controls > threshold)
    }
    else if (direction == "<") {
        tp <- sum(cases >= threshold)
        tn <- sum(controls < threshold)
    }
    
    return(c(sp = tn/length(controls), se = tp/length(cases)))
}



perfs.fun <- function(x, thresholds, condition) {
    
    if (condition == "0") {
        res <- as.numeric(x[1]>=thresholds[2] | (x[1]>= thresholds[1] & x[2] >= thresholds[3]) )
    } else if (condition == "1") {
        res <- as.numeric(x[1]<thresholds[1] | (x[1] < thresholds[2] & x[2] < thresholds[3]) )
    }
    
    return(res)
}



perfs.speed=function (thresholds, controls, cases, direction) 
{
    
    # calculate the sensitivity and specificity
    if (direction == "<") {
        tp <- apply(cases, 1, perfs.fun, thresholds=thresholds, condition="0")
        tn <- apply(controls, 1, perfs.fun, thresholds=thresholds, condition="1")
    }
    else if (direction == ">") {
        tp <- apply(cases, 1, perfs.fun, thresholds=thresholds, condition="1")
        tn <- apply(controls, 1, perfs.fun, thresholds=thresholds, condition="0")
    }
    
    sp=rowMeans(tn)
    se=rowMeans(tp)
    rm(tp);rm(tn)
    # gc()
    
    return(data.frame(sp=sp, se=se))
}


upperFinder <- function(ss){
    
    ## type: upper
    uss=unique(ss)
    
    uss.upper1=aggregate(uss[,2], by = list(uss[,1]), max)
    uss.upper2=aggregate(uss[,1], by = list(uss[,2]), max)
    uss.upper2 = uss.upper2[,c(2,1)]
    
    colnames(uss.upper1) = c("sp", "se")
    colnames(uss.upper2) = c("sp", "se")
    
    uss.Upper = unique(rbind(uss.upper1,uss.upper2))
    uss.Upper = uss.Upper[order(uss.Upper[,2],decreasing=T),]
    
    ss.upper = as.matrix(uss.Upper)
    
    for (i in 1:(nrow(uss.Upper)-1)) {
        tbl22 = uss.Upper[c(i,i+1),]
        
        if (!all(tbl22[1,1] == tbl22[,1]) & !all(tbl22[1,2] == tbl22[,2])) {
            ss.upper = rbind(ss.upper,c(min(tbl22[,1]), min(tbl22[,2])))
        } 
    }
    
    ss.upper = ss.upper[order(ss.upper[,1],decreasing=F),]
    ss.upper = ss.upper[order(ss.upper[,2],decreasing=T),]

    return(ss.upper)
}


auc.nocheck<- function(ss) {
    
    # if (is.unsorted(ss[,1])) {
    #     ss[,1] <- rev(ss[,1])
    #     ss[,2] <- rev(ss[,2])
    # }
    
    sp = ss[,1]
    se = ss[,2]
    diffs.x <- sp[-1] - sp[-length(sp)]
    means.vert <- (se[-1] + se[-length(se)])/2
    auc <- sum(means.vert * diffs.x)
    return(auc)
}


# wrapup function
auc.speed <- function(controls, cases, direction = c("auto", "<",">")) {
    
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

    kk = perfs.speed(thresholds=tt,controls=controls, cases=cases, direction=direction)
    
    # obtain upper curve
    ss.upper = upperFinder(kk)
    
    # calculate auc
    auc.upper = auc.nocheck(ss.upper)
    
    return(auc.upper)
}


roc.new <- function(controls, cases, direction = c("auto", "<",">"), 
                    auc=T) 
{
    
    # set direction by median
    direction <- match.arg(direction)
    
    if (direction == "auto" && median(controls[,1]) <= median(cases[,1])) {
        direction <- "<"
    } else if (direction == "auto" && median(controls[,1]) > median(cases[,1])) {
        direction <- ">"
    } 
    
    # create legit threshold candidate
    tmat=thresholds.new(controls,cases)
    tt=expand.grid(tmat[,1],tmat[,2],tmat[,3])
    tt=tt[tt[,1]<=tt[,2],]
    
    ### single step
    st=tmat[,1]
    
    ss.single= t(sapply(st, perfs.single, controls = controls, cases = cases, direction = direction))
    
    
    # sensitivity and specificity
    kk = perfs.speed(thresholds=tt,controls=controls, cases=cases, direction=direction)
    
    # constructing results object
    fnltbl = cbind(tt,kk)
    
    colnames(fnltbl)[c(1:3)] <- c("theta1", "theta2", "theta3")
    
    # different type of curves (sensitivity and specificity)
    ss.upper = upperFinder(kk[,c(1:2)])
    
    
    # create roca object
    roc <- list()
    class(roc) <- "mtroc"
    
    roc$ss.raw <- kk
    
    roc$ss.single <- ss.single
    
    roc$direction <- direction
    roc$cases <- cases
    roc$controls <- controls
    
    roc$results <- fnltbl
    
    roc$ss.upper <- ss.upper
    
    if (auc) 
        roc$auc <- auc.nocheck(ss.upper)
    
    return(roc)
}

# roc <- roc.new(controls=controls, cases=cases)
