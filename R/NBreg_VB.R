NBreg_VB <- function(counts, X, cond, maxIter = 1000)
{
    y <- as.matrix(counts)
    dd <- dim(y)[1]
    idx.nz <- rowSums(y)!=0
    y <- y[idx.nz,]
    K <- dim(y)[1]       
    J <- dim(y)[2]     
    P <- dim(X)[1]
    idx.group1 <- which(cond==1)
    idx.group2 <- which(cond==2)
    # hyperparameters
    a0 <- b0 <- c0 <- d0 <- e0 <- f0 <- g0 <- 0.01
    
    # variational parameters
    rj <- rep(1,J)             # \tilde{r}_j
    aj <- hj <- rep(1,J)
    cp <- dp <- rep(1,P)
    cp <- (c0+K/2)*cp
    muk <- matrix(0,P,K)
    sigmak <- vector("list",K)
    b <- a0 + b0
    g <- 1
    
    # mean quantities
    Elogr <- digamma(aj)-log(hj)
    EL <- (rep(1,K) %*% t(rj)) * (digamma(y+rep(1,K) %*% t(rj))- rep(1,K) %*% t(digamma(rj)))
    Er <- aj/hj
    Ealpha <- cp/dp
    Eh <- b/g
    Elog1pexp <- Ew <- matrix(0,K,J)
    for(k in 1:K){
        betak <- rnorm(P)
        Elog1pexp[k,] <- logOnePlusExp(betak %*% X)
        Ew[k,] <- (y[k,]+Er)*(tanh(betak %*% X/2)/(2*betak %*% X))
    }
    
    for(iter in 1:maxIter)
    {
        # update L_{jk}
        rj <- exp(Elogr)
        # r_{j}
        aj <- a0 + colSums(EL)
        hj <- Eh + colSums(Elog1pexp)
        Er <- aj/hj
        Elogr <- digamma(aj)-log(hj)
        # \beta_k
        temp <- 0
        for(k in 1:K){
            sigmak[[k]] <- solve(X %*% diag(Ew[k,]) %*% t(X)+diag(Ealpha), diag(P))
            muk[,k] <- sigmak[[k]] %*% (0.5*X %*% (y[k,]-Er))
            # update mean quantities
            betak <- mvrnorm(n = 500, mu = muk[,k], Sigma = sigmak[[k]], tol = 1e-30)    # 300 by P
            Elog1pexp[k,] <- colMeans(logOnePlusExp(betak %*% X))
            Ew[k,] <- (y[k,]+Er)*colMeans(tanh(betak %*% X/2)/(2*betak %*% X))
            temp <- temp + diag(sigmak[[k]]) + muk[,k]^2
        }
        # \alpha_p
        for (k in 1)
        dp <- d0 + 0.5*temp
        Ealpha <- cp/dp
        # h
        g <- g0 + sum(Er)
        Eh <- b/g
    }
    #mu <- matrix(0,P,dd)
    #mu[,idx.nz] <- muk
    return(list(mu=muk,sigma=sigmak,r=Er,h=Eh,alpha=Ealpha))
}