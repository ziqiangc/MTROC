# FUN empirical FPF/TPF

source("FUN_GreyZone_v2.R")


sim.fun <- function(d, var, rho, n){
    Sigma <- matrix(c(var,rho,rho,var),2,2)
    controls=mvrnorm(n, mu=c(0,0), Sigma = Sigma)
    cases=mvrnorm(n, mu=c(d,d), Sigma = Sigma)
    
    return(list(controls=controls, cases=cases))
}


empirical.sim <- function(d, var, rho, n, theta){
    
    simdat <- sim.fun(d=d, var=var, rho=rho, n=n)
    # theta.mat <- theta.gen(d=d, var=var, k=2, length=length)
    
    pfs <- perfs.new(controls=simdat$controls, cases=simdat$cases, direction = "<", thresholds=theta)
    
    
    return(pfs)
}

