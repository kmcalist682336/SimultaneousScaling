#Sample the top level latent scores
#Relatively simple
#Not a lot of speed to gain, here
#Can all be evaluated in parallel, but the cost doesn't really outweigh the gain
#Arguments: yy is matrix (p x n) of augmented vote data
#Lambda is the p x K matrix of top level factor loadings
#Alpha is the p-vector of item level intercepts
#n is the number of voters
#K is the number of latent dimensions at the current step
sample.omega <- function(yy,lambda,alpha,n,K){
  v.d <- solve((t(lambda)%*%lambda) + diag(1,K))
  v.dd <- v.d%*%t(lambda)
  draw.omg <- function(i){
    ya <- yy[,i] + alpha
    v.ddd <- v.dd%*%ya
    n.ddd <- rmvnorm(1,mean = c(v.ddd),sigma = v.d)
    return(as.matrix(n.ddd))
  }
  draw.omg <- Vectorize(draw.omg,"i")
  new.omg <- draw.omg(i = 1:n)
  return(t(new.omg))
}