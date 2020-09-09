# FUN parametric FPF/TPF

sim.fun <- function(d, var, rho, n){
    Sigma <- matrix(c(var,rho,rho,var),2,2)
    controls=mvrnorm(n, mu=c(0,0), Sigma = Sigma)
    cases=mvrnorm(n, mu=c(d,d), Sigma = Sigma)
    mu.c <- apply(controls, 2, mean)
    mu.d <- apply(cases, 2, mean)
    sigma.c <- cov(controls)
    sigma.d <- cov(cases)
    
    return(list(mu.c=mu.c, mu.d=mu.d, sigma.c=sigma.c, sigma.d=sigma.d))
}


pf.fun <- function(mu, sigma, theta){
    p1=pnorm(q = (mu[1]-theta[2])/sqrt(sigma[1,1]) )
    p2=pmvnorm(lower = c(theta[1], theta[3]), upper = c(theta[2], Inf), mean = mu, sigma = sigma)
    p=p1+p2
    return(p)
}



param.sim <- function(d, var, rho, n, theta){
    
    params <- sim.fun(d=d, var=var, rho=rho, n=n)
    # theta.mat <- theta.gen(d=d, var=var, k=2, length=length)
    
    # fpf=apply(theta.mat, 1, pf.fun, mu=params$mu.c, sigma=params$sigma.c)
    # tpf=apply(theta.mat, 1, pf.fun, mu=params$mu.d, sigma=params$sigma.d)
    
    fpf=pf.fun(mu=params$mu.c, sigma=params$sigma.c, theta=theta)
    tpf=pf.fun(mu=params$mu.d, sigma=params$sigma.d, theta=theta)
    
    pfs <- c(fpf,tpf)
    
    
    return(pfs)
}
