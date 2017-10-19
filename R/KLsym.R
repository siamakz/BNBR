KLsym <- function(x1, x2)
{
    n <- nrow(x1)
    if (nrow(x2)!=n) stop("Number of rows should be the same.")
    nbins <- 100
    c <- 1.5
    KLdist <- NULL
    for (i in 1:n)
    {
        #find outliers
        q <- quantile(c(x1[i, ], x2[i, ]))
        Q1 <- q[2]
        Q3 <- q[4]
        IQD <- Q3 - Q1
        xmin <- max(Q1 - c*IQD,0)
        xmax <- Q3 + c*IQD
        idx <- which(x1[i,]>=xmin & x1[i,]<=xmax)
        f1 <- (hist(x1[i, idx], breaks = seq(xmin, xmax, length.out = nbins), plot = FALSE)$counts)/length(idx) + 1e-10
        idx <- which(x2[i,]>=xmin & x2[i,]<=xmax)
        f2 <- (hist(x2[i, idx], breaks = seq(xmin, xmax, length.out = nbins), plot = FALSE)$counts)/length(idx) + 1e-10
        KLdist <- c(KLdist, sum((f1-f2)*log(f1/f2)))
    }
    return(KLdist)
}