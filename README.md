# Bayesian Negative Binomial Regression for Differential Expression Analysis

This is a R package for differential expression analysis of seqencing count data. The main function NBregDE2 takes as input:

*  **counts**:  a matrix of counts, rows are corresponding to genes and columns are corresponding to samples
*  **X**: design matrix
*  **cond**: the index of covariate corresponding to the main treatment
*  **idx.cond**: a vector of covariate indices that the effect of their combination on gene expression is under study
*  **Burnin**: number of burn-in iterations in MCMC
*  **Collections**: number of collected posterior samples after burn-in
*  **PGTruncation**: the truncation level used for genrating random numbers from Polya-Gamma distribution
*  **randtry**: to be used for set.seed()

## Example 

Given the gene count matrix and two covariates as _sex_ and _cond_, the ranking of genes with respect to _cond_ can be obtained as:

```R
Burnin <- 1000
Collections <- 1000
X <- as.matrix(t(model.matrix(~sex+cond)))
idx.cond <- 3
res <- NBregDE2(dat,X,cond,idx.cond, Burnin,Collections)
res$kl            # the KL-divergence score of genes, high values correspond to high DE
res$beta          # posterior means of regression coefficient for the covariate of interest (i.e., cond)
```
