# I^2 for Multilevel and Multivariate Models
"Now we will consider the same type of generalization, but for a multivariate model 
with non-independent sampling errors. Therefore, not only do we need to account for
heterogeneity and dependency in the underlying true effects, but we also now need to
specify covariances between the sampling errors. Two things are worth noting here. 
First of all, we allow the amount of heterogeneity to differ for all the outcomes
by using an unstructured variance-covariance matrix for the true effects (i.e., struct='UN').
Second, V is the variance-covariance matrix of the sampling errors, which is no longer 
diagonal. The W  matrix described earlier is actually the inverse of the variance-covariance
matrix of the sampling errors.
This statistic can be thought of as the overall I2 value that indicates how much of the total 
variance can be attributed to the total amount of heterogeneity."

calculate_I2 <- function(y, V){
  W <- solve(V)
  X <- model.matrix(y)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2 <- 100 * y$tau2 / (y$tau2 + (y$k-y$p)/sum(diag(P)))
  
  return(I2)
}