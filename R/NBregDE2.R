#' Bayesian Negative Binomial regression for differential expression analysis of
#' sequencing count data with confounding factors
#'
#' @param counts a matrix of counts, rows are corresponding to genes and columns
#' are corresponding to samples.
#' @param X design matrix
#' @param cond the index of covariate corresponding to the main treatment
#' @param idx.cond a vector of covariate indices that the effect of their combination
#' on gene expression is under study
#' @param Burnin Number of burn-in iterations in MCMC
#' @param Collections Number of collected posterior samples after burn-in
#' @param PGTruncation the truncation level used for genrating random numbers
#' from Polya-Gamma distribution
#' @param randtry to be used for set.seed()
#' @export
NBregDE2 <- function(counts, X, cond, idx.cond, Burnin = 1000L, Collections = 1000L, PGTruncation = 60L, randtry = 2017)
{
    set.seed(randtry)
    y <- as.matrix(counts)
    idx.nz <- rowSums(y)!=0
    y <- y[idx.nz,]
    ngenes <- dim(y)[1]       # K
    nsamples <- dim(y)[2]     # J
    idx.group1 <- which(cond==1)
    idx.group2 <- which(cond==2)
    # hyperparameters
    a0 <- b0 <- c0 <- d0 <- e0 <- f0 <- g0 <- 0.01

    # coefficients matrix beta, P by K
    P <- dim(X)[1]
    Beta <- matrix(0, P, ngenes)
    r <- rep(100, nsamples)
    h <- 1
    alpha <- rep(1, P)
    ####
    ter <- NBreg_VB(counts,X,cond,100)
    Beta <- ter$mu
    r <- ter$r
    alpha <- ter$alpha
    h <- ter$h
    ###

    Psi <- matrix(0, ngenes, nsamples)
    XYT <- X %*% t(y)

    beta.means <- matrix(0, P, ngenes)

    beta2.samples <- beta3.samples <- matrix(0, ngenes, Collections)
    theta.means <- matrix(0, ngenes, 2*Collections)
    iterMax <- Burnin+Collections
    yy <- pmin(y,10000)

    for (iter in 1:iterMax)
    {
        cat(iter, '\n')

        # Sample r_j
        ell <- rep(0, nsamples)
        for (j in 1:nsamples)
        {
            yr <- y[y[,j]>10000,j]
            ell[j] <- CRT_sum(yy[,j], r[j]) + rpois(1, r[j]*(sum(digamma(yr+r[j]))-length(yr)*digamma(10000+r[j])))
        }
        r <- rgamma(nsamples, a0 + ell, rate = h + colSums(logOnePlusExp(Psi)))

        # Sample alpha
        alpha <- rgamma(P, c0 + ngenes/2, rate = d0 + 0.5*rowSums(Beta^2))

        # Smple h
        h <- rgamma(1, b0+nsamples*a0, rate = g0+sum(r))


        # Sample omega
        temp <- y+matrix(rep(r, ngenes), ngenes, nsamples, byrow = T)
        omega <- PolyaGamRnd_Gam(c(temp), c(Psi), Truncation = PGTruncation)
        omega <- matrix(omega, ngenes, nsamples)  # K by J

        # Sample Beta, phi and Psi
        A <- diag(alpha)
        for (k in 1:ngenes)
        {
            if (any(eigen(A + X %*% diag(omega[k,]) %*% t(X))$values<=0)) stop('ter')
            temp <- solve(chol(A + X %*% diag(omega[k,]) %*% t(X)))
            Beta[,k] <- temp %*% (rnorm(P) + t(temp) %*% (0.5*(XYT[,k]-X %*% r)))
        }

        if (iter>Burnin)
        {
             beta.means <- beta.means + Beta
#            beta2.samples[,iter-Burnin] <- Beta[2,]
#            beta3.samples[,iter-Burnin] <- Beta[3,]
#             temp <- 1 + exp(t(Beta) %*% X)
             theta.means[, iter-Burnin] <- exp(Beta[1,])
             if (length(idx.cond)==1){
                theta.means[, iter-Burnin+Collections] <- exp(Beta[1,]+Beta[idx.cond,])
             }
             if (length(idx.cond)>1){
                theta.means[, iter-Burnin+Collections] <- exp(Beta[1,]+colSums(Beta[idx.cond,]))
             }
        }
    }

if (0){
    out1 <- out2 <- matrix(0, dim(counts)[1], Collections)
    out1[idx.nz, ] <- beta2.samples
    out2[idx.nz, ] <- beta3.samples
    out <- matrix(0, dim(counts)[1], 2*Collections)
    out[idx.nz, ] <- theta.means
   }
   kl <- KLsym(theta.means[,1:Collections],theta.means[,(Collections+1):(2*Collections)])
   out1 <- rep(0,dim(counts)[1])
   out1[idx.nz] <- kl
   out <- matrix(0, P, dim(counts)[1])
   out[,idx.nz] <- beta.means/Collections
#    return(list(beta2 = out1, beta3 = out2, theta = out))
    return(list(beta=out,kl=out1))
}
