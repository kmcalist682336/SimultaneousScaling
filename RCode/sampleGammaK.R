#Sample the factor loading precisions
#TOp level loading precisions are assumed to be common over all p values in the kth dim
#Doesn't make much difference to assume different variances
#Simple, no changes in future
#Arguments: lambda is matrix (p x K) of factor loadings
#K is the number of current dimensions on the top level
sample.gamma.k <- function(lambda,K){
  mk <- apply(sign(abs(lambda)),2,sum)
  sls <- apply(lambda^2,2,sum)
  new.gamma.k <- c()
  for(i in 1:K){
    new.gamma.k[i] <- rgamma(1,mk[i],sls[i])
  }
  return(new.gamma.k)
}