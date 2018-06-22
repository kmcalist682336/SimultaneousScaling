sample.zz.lam <- function(yy,lambda,zz,alpha,omega,pis,gamma.k,n,p,K){
  #Sample the matrices zz and lambda
  #Return list with two objects: zz and lambda
  #Everything should be calculated in sequence
  #1-p over k then 1-K
  #Fast computation, no tricks needed
  new.lambda <- lambda
  new.zz <- zz
  #For loops is really the best way to approach this
  #Julia this function later
  #First, precompute w'w + y
  gam <- c()
  for(i in 1:K){
    gam[i] <- (t(as.matrix(omega[,i]))%*%as.matrix(omega[,i])) + gamma.k[i]
  }
  for(j in 1:K){
    for(i in 1:p){
      lam.jk <- as.matrix(new.lambda[i,])
      lam.jk[j,1] <- 0
      alpha.j <- as.matrix(rep(alpha[i,1],n)) 
      ehat <- as.matrix(yy[i,]) - t((t(lam.jk)%*%t(omega))) + alpha.j
      mu <- (1/gam[j])*(t(as.matrix(omega[,j]))%*%ehat)
      l1 <- (.5*log(gamma.k[j])) - (.5*log(gam[j])) + (.5*gam[j]*mu)
      l2 <- log(pis[i,j]) - log(1 - pis[i,j])
      ll <- l1 + l2
      if(ll > 10){
        nz <- 1
      }else{
        if(ll < -10){
          nz <- 0
        }else{
          pv <- exp(ll)/(1 + exp(ll))
          nz <- rbinom(1,1,pv)
        }
      }
      if(nz == 0){
        nl <- 0
      }else{
        nl <- rnorm(1,mu,sqrt(1/gam[j]))
      }
      new.lambda[i,j] <- nl
      new.zz[i,j] <- nz
    }
  }
  out.list <- list()
  out.list[[1]] <- new.zz
  out.list[[2]] <- new.lambda
  names(out.list) <- c("zz","lambda")
  return(out.list)
}