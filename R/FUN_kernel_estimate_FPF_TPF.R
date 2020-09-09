# # FUN kernel FPF/TPF


source("FUN_kcde_wrapup.R")

sim.fun <- function(d, var, rho, n){
    Sigma <- matrix(c(var,rho,rho,var),2,2)
    controls=mvrnorm(n, mu=c(0,0), Sigma = Sigma)
    cases=mvrnorm(n, mu=c(d,d), Sigma = Sigma)
    
    return(list(controls=controls, cases=cases))
}



pf.fun <- function(Fhat, theta) {
    p1=kcde.wrapup(Fhat, eval.points=c(theta[2]))$estimate
    m2=matrix(c(theta[1],theta[3],theta[2],theta[3]),nrow=2,byrow = T)
    p2=kcde.wrapup(Fhat, eval.points=m2)$estimate
    p=p1+(p2[1]-p2[2])
    return(p)
}


nparam.sim <- function(d, var, rho, n, theta){
    
    simdat <- sim.fun(d=d, var=var, rho=rho, n=n)
    # theta.mat <- theta.gen(d=d, var=var, k=2, length=length)
    
    Fhat1 <- kde.wrapup(x=simdat$controls, supp = 3.7)
    Fhat2 <- kde.wrapup(x=simdat$cases, supp = 3.7)
    
    # Fhat1 <- kde.wrapup(x=simdat$controls, supp = max(theta.mat))
    # Fhat2 <- kde.wrapup(x=simdat$cases, supp = 2*max(theta.mat))
    
    # fpf=mclapply(as.data.frame(t(theta.mat)), pf.fun, Fhat=Fhat1, mc.cores=4)
    # tpf=mclapply(as.data.frame(t(theta.mat)), pf.fun, Fhat=Fhat2, mc.cores=4)
    
    fpf=pf.fun(Fhat=Fhat1, theta=theta)
    tpf=pf.fun(Fhat=Fhat2, theta=theta)
    
    
    
    pfs <- c(fpf,tpf)
    
    return(pfs)
}
