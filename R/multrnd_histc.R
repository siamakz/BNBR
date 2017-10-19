multrnd_histc <- function(n,p)
{
    edges <- c(0, cumsum(p))
    l <- length(edges)
    r <- hist(runif(n)*edges[l], edges, plot = FALSE)$counts
    return(r)
}