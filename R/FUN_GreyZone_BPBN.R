###
# change perf.single
#
###



# function to find cutoff points (2-step)
# can be modified to retrun a list (# of theta2 is not necessary same as theta1)

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
    return(list(theta1.lower,theta1.upper,theta2))
}



# function to calculate the sensitivity and specificity

perfs.new=function (thresholds, controls, cases, direction) 
{
    # calculate the sensitivity and specificity
    if (direction == "<") {
        tp <- sum(  ifelse(cases[,1] < thresholds[1], 0, 
                           ifelse(cases[,1] >= thresholds[2], 1, 
                                  ifelse(cases[,2] < thresholds[3], 0, 1))) )
        tn <- sum(  ifelse(controls[,1] < thresholds[1], 1, 
                           ifelse(controls[,1] >= thresholds[2], 0, 
                                  ifelse(controls[,2] < thresholds[3], 1, 0))) )
    }
    else if (direction == ">") {
        tp <- sum(  ifelse(cases[,1] < thresholds[1], 1, 
                           ifelse(cases[,1] >= thresholds[2], 0, 
                                  ifelse(cases[,2] < thresholds[3], 1, 0))) )
        tn <- sum(  ifelse(controls[,1] < thresholds[1], 0, 
                           ifelse(controls[,1] >= thresholds[2], 1, 
                                  ifelse(controls[,2] < thresholds[3], 0, 1))) )
    }
    
    
    # calculate the probability of requiring the 2nd step
    
    # p <- mean(c(controls[,1], cases[,1]) >= thresholds[1] & c(controls[,1], cases[,1]) < thresholds[2])
    
    return(c(sp = tn/length(controls[,1]), se = tp/length(cases[,1])))
}


## single step
perfs.single=function (threshold, controls, cases, direction) 
{
    # controls = controls[,1]
    # cases = cases[,1]
    
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



# obtain different type curves
curveFinder <- function(ss){
    
    ## type: unique
    # ss.unique = unique(ss)
    
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
    
    
    ## type: lower
    uss.lower1=aggregate(uss[,2], by = list(uss[,1]), min)
    uss.lower2=aggregate(uss[,1], by = list(uss[,2]), min)
    uss.lower2 = uss.lower2[,c(2,1)]
    
    colnames(uss.lower1) = c("sp", "se")
    colnames(uss.lower2) = c("sp", "se")
    
    uss.Lower = unique(rbind(uss.lower1,uss.lower2))
    uss.Lower = uss.Lower[order(uss.Lower[,2],decreasing=T),]
    
    ss.lower = as.matrix(uss.Lower)
    
    for (i in 1:(nrow(uss.Lower)-1)) {
        tbl22 = uss.Lower[c(i,i+1),]
        
        if (!all(tbl22[1,1] == tbl22[,1]) & !all(tbl22[1,2] == tbl22[,2])) {
            ss.lower = rbind(ss.lower,c(max(tbl22[,1]), max(tbl22[,2])))
        } 
    }
    
    ss.lower = ss.lower[order(ss.lower[,2],decreasing=T),]
    
    
    ss.list <- list(ss.unique=uss,
                    ss.upper=ss.upper,
                    ss.lower=ss.lower)
    
    return(ss.list)
}



# function to obtain auc
auc.roca <- function(roca, type.auc=c("single", "sort", "upper", "lower", "average")) 
{
    if(type.auc=="single") ss = roca$ss.single
    if(type.auc=="sort") ss = roca$ss.list$ss.sort
    if(type.auc=="upper") ss = roca$ss.list$ss.upper
    if(type.auc=="lower") ss = roca$ss.list$ss.lower
    
    
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




##########test##############

roc.new <- function(controls, cases, direction = c("auto", "<",">"), strategy = c("BE", "BN", "BP"),
                    auc=T, type.auc=c("single","sort", "upper", "lower")) 
{
    
    # set direction by median
    direction <- match.arg(direction)
    
    if (direction == "auto" && median(controls[,1]) <= median(cases[,1])) {
        direction <- "<"
    } else if (direction == "auto" && median(controls[,1]) > median(cases[,1])) {
        direction <- ">"
    } 
    
    # create legit threshold candidate
    tlst=thresholds.new(controls,cases)
    if (strategy =="BE") {
        tt=expand.grid(tlst[[1]],tlst[[2]],tlst[[3]])
        tt=tt[tt[,1]<=tt[,2],]
    } else if (strategy == "BN") {
        tt=expand.grid(tlst[[1]],Inf,tlst[[3]])
    } else if (strategy == "BP") {
        tt=expand.grid(-Inf,tlst[[2]],tlst[[3]])
    }


    
    ### single step
    st1 = tlst[[1]]
    
    st2 = tlst[[3]]
    
    ss.single1= t(sapply(st1, perfs.single, controls = controls[,1], cases = cases[,1], direction = direction))
    ss.single2= t(sapply(st2, perfs.single, controls = controls[,2], cases = cases[,2], direction = direction))
    
    # sensitivity and specificity
    kk=apply(tt, 1, perfs.new, controls=controls, cases=cases, direction=direction)
    
    # constructing results object
    fnltbl = cbind(tt,t(kk))
    
    colnames(fnltbl)[c(1:3)] <- c("theta1", "theta2", "theta3")
    
    # different type of curves (sensitivity and specificity)
    ss.list = curveFinder(t(kk)[,c(1:2)])
    
    
    # create roca object
    roc <- list()
    class(roc) <- "roca"
    
    roc$ss.raw <- t(kk)[,c(1:2)]
    
    roc$ss.single1 <- ss.single1
    roc$ss.single2 <- ss.single2
    
    roc$direction <- direction
    roc$cases <- cases
    roc$controls <- controls
    
    roc$results <- fnltbl
    
    roc$ss.list <- ss.list
    
    # if (auc) 
    #     roc$auc <- auc.roca(roc, type.auc=type.auc)
    
    return(roc)
}


#### plot function ####

plot.roca <- function(roca, auc=T)
{
    
    ss.single = roca$ss.single
    ss.raw = roca$ss.raw
    ss.upper = roca$ss.list$ss.upper
    ss.lower = roca$ss.list$ss.lower
    
    smoothScatter(1-ss.raw[,1], ss.raw[,2], pch = 21, nrpoints=length(ss.raw[,1]), 
                  colramp = colorRampPalette(c("white", blues9[1:8])), 
                  xlim=c(0,1), ylim=c(0,1), xlab = "FPF", ylab = "TPF", cex=0.6)
    lines(1-ss.single[,1], ss.single[,2], lty=1, lwd=3)
    lines(1-ss.upper[,1], ss.upper[,2], col="darkred", lty=2, lwd=2)
    lines(1-ss.lower[,1], ss.lower[,2], col="darkgreen", lty=2, lwd=2)
    
    if (auc) {
        legend("bottomright", legend=c(paste("Regular ROC", "(AUC=",auc.roca(roc,type.auc=c("single")),")"), 
                                       paste("Greyzone Upper", "(AUC=",auc.roca(roc,type.auc=c("upper")),")"),
                                       paste("Greyzone Lower", "(AUC=",auc.roca(roc,type.auc=c("lower")),")")), 
               lty=c(1,2,2,2), 
               col=c("black", "darkred", "darkgreen"), lwd=2, bty="n")
        
    } else {
        legend("bottomright", legend=c("Regular ROC", "Greyzone Upper", "Greyzone Lower"), lty=c(1,2,2), 
               col=c("black", "darkred", "darkgreen"), lwd=2, bty="n")
        
    }
    
}






# roc=roc.new(controls, cases, direction="<", auc=T, type.auc=c("single"))
# 
# 
# png("Figures/roc_comparison_test.png", width = 6, height = 6, units = "in", res = 300)
# plot(roc)
# dev.off()





