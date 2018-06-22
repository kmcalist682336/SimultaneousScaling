#Sample item level intercepts for the top level
#Intercepts are meaningful here because there is theory as to what the intercepts mean
#Recall the top level model: y_i,j = lambda_j omega_i - alpha_j
#Simple evaluation, no improvements needed
#Maybe allow prior specification here later, but its hard to imagine that anyone has meaningful priors on this parameter
#Assume flat, improper normal prior
#Arguments: yy is augmented vote data (p x n)
#Lambda is top level factor loadings (p x K)
#Omega is top level factor scores (n x K)
#p is number of items/bills/documents
sample.alpha <- function(yy,lambda,omega,p){
  pred.mat <- lambda%*%t(omega)
  pred.mat <- pred.mat - yy
  pp.mat <- (1/n)*apply(pred.mat,1,sum)
  new.alpha <- c()
  for(i in 1:p){
    new.alpha[i] <- rnorm(1,pp.mat[i],sqrt(1/n))
  }
  return(as.matrix(new.alpha))
}