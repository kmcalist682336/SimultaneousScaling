###Create Sim Data for Multilevel RollCall model
###Create Utility functions for MCMC
rm(list=ls())
set.seed(12345)
library(mvtnorm)
library(truncnorm)
library(logisticPCA)
library(abind)
library(doParallel)
library(foreach)
library(data.table)
library(generalizedPCA)
library(Brobdingnag)
library(parallelDist)
library(amap)
library(parallel)
library(rootSolve)
euc_dist <- function(m) {mtm <- Matrix::tcrossprod(m); sq <- rowSums(m^2);  sqq <- outer(sq,sq,"+") - (2*mtm); sqq[sqq <0] <- 0; return(sqq)}
#####Generate true values
###Start with bottom level
###Use bottom level to determine top level
#n is the number of observations
#p is the number of bills
#d is the number of tokens
#l is the number of text dims
#k is the number of vote dimes
n <- 500
p <- 100
d <- 100
L <- 3
K <- 3
#gen eta.l
t.eta.l <- as.matrix(runif(L,.1,.3))
#gen R
t.r <- matrix(ncol = L, nrow = p)
for(i in 1:p){
  for(j in 1:L){
    t.r[i,j] <- rbinom(1,1,t.eta.l[j])
  }
}
#gen A
t.a <- matrix(ncol = L, nrow = p)
for(i in 1:p){
  for(j in 1:L){
    t.a[i,j] <- sample(c(-1,1),1)*runif(1,.1,2)
  }
}
t.a <- t.a*t.r
#gen B
t.b <- matrix(ncol = L, nrow = d)
for(i in 1:d){
  for(j in 1:L){
    t.b[i,j] <- runif(1,-2,2)
  }
}

#gen augmented dats, M
mm <- t.b%*%t(t.r*t.a)
#Use mm to generate observed mstar - poisson vars
mstar <- mm
for(i in 1:d){
  mms <- cbind(seq(1,p),mstar[i,])
  mms <- mms[order(mms[,2]),]
  ii <- rpois(1,10)
  word.vals <- sort(rpois(p,ii))
  mms <- cbind(mms,word.vals)
  mms <- mms[order(mms[,1]),]
  mstar[i,] <- mms[,3]
}
#generate the true dist mat for a
#abs of location for now
make.norm.dist.mat <- function(m){
  t.dist <- euc_dist(abs(m))
  t.dist <- exp(-as.matrix(t.dist))
  sms <- apply(t.dist,1,sum)
  sms <- matrix(ncol = p,rep(sms,p))
  t.dist <- t.dist/sms
  return(t.dist)
}
t.dist <- make.norm.dist.mat(t.a)
#Map each of the l dims to the k and make pi.star that way
l.to.k <- sample(seq(1,K),L,replace = T)
t.pi.star <- matrix(ncol = K, nrow = p,runif(K*p,0,.05))
for(i in 1:p){
  rr <- sort(unique(t.r[i,]*l.to.k))[-1]
  for(j in 1:length(rr)){
    t.pi.star[i,rr[j]] <- runif(1,.8,1)
  }
}
t.pi <- matrix(ncol = K, nrow = p)
for(i in 1:p){
  for(j in 1:K){
    t.pi[i,j] <- sum(t.dist[i,]*t.pi.star[,j])
  }
}
t.z <- matrix(ncol = K, nrow = p)
for(i in 1:p){
  for(j in 1:K){
    t.z[i,j] <- rbinom(1,1,t.pi[i,j])
  }
}
t.lam <- matrix(ncol = K, nrow = p)
for(i in 1:p){
  for(j in 1:K){
    t.lam[i,j] <- sample(c(-1,1),1)*runif(1,.1,2)
  }
}
t.lam <- t.lam*t.z
t.omg <- matrix(ncol = K, nrow = n,rnorm(K*n))
t.alp <- as.matrix(rnorm(p,0,1))
t.yy <- (t.lam%*%t(t.omg)) - matrix(rep(t.alp,n),ncol = n)
ystar <- sign(sign(t.yy) + 1)
#Randomly insert missingness in the ystar matrix to callibrate functions for missingness
#Note that there can be missingness in the vote matrix, but not the text matrix
na.mat <- matrix(sample(c(0,1),n*p,replace = T, prob = c(.9,.1)), ncol = n, nrow = p)
for(i in 1:p){
  for(j in 1:n){
    if(na.mat[i,j] == 1){
      ystar[i,j] <- NA
    }
  }
}
########################################################################################################
#ystar is the binary vote data
#mstar is the count data for the text tokens
#Create initial values and containers for quantities of interest
#Values to keep:
#alpha - p
#omega - K x n
#lambda - p x K
#z - p x K
#pi - p x K
#pi* - p x K
#eta - K
#gamma - K
#A - L x p
#Kappa - p x p
#R - L x p
#xi - L
#B - d x L
#sigma - d x L
#Set the initial values for K+ and L+
L <- 25
K <- 25
#Let's set n,p, and d from the observed values
n <- dim(ystar)[2]
p <- dim(ystar)[1]
d <- dim(mstar)[1]
#Randomly generate starting values for omega, lambda, and alpha
#Set all z=1 to start
omega <- matrix(ncol = K, nrow = n, rnorm(n*K,0,1))
lambda <- matrix(ncol = K, nrow = p, rnorm(p*K,0,1))
zz <- matrix(ncol = K, nrow = p, c(1))
alpha <- as.matrix(rnorm(p,0,.01))
#Generate starting values for pi.star
pi.star <- matrix(ncol = K, nrow = p, rbeta(K*p,1,1))
#Generate k starting values for gamma and eta
gamma.k <- as.matrix(runif(K,.1,2))
eta.k <- as.matrix(rbeta(K,1,1))
#Write utility function to return augmented set of data given ystar
augment.vote.data <- function(omega,lambda,alpha,ystar,n,p){
  #Start by matrixing alpha
  #P items repeated N times
  #alpha rbind n times
  alpha.mat <- matrix(rep(alpha,n),ncol = n)
  pred.mat <- data.table((lambda%*%t(omega)) - alpha.mat)
  setnames(pred.mat,old = names(pred.mat), new = as.character(seq(1,n)))
  pred.mat <- pred.mat[,id := seq(1,p)]
  ys <- data.table(ystar)
  setnames(ys,old = names(ys), new = as.character(seq(1,n)))
  ys <- ys[,id := seq(1,p)]
  pred.mat <- reshape2::melt(pred.mat, id.var = "id")
  ys <- reshape2::melt(ys, id.var = "id")
  setkeyv(x = pred.mat,cols = c("id","variable"))
  setkeyv(x = ys,cols = c("id","variable"))
  js <- pred.mat[ys, nomatch=0]
  rm(ys,pred.mat)
  setnames(js, names(js), c("id","vars","pred.m","obs.v"))
  js.list <- list()
  js.list[[1]] <- js[is.na(obs.v) == T]
  js.list[[2]] <- js[obs.v == 0]
  js.list[[3]] <- js[obs.v == 1]              
  rm(js)
  gc()
  js.list[[1]] <- js.list[[1]][,new.aug := rnorm(dim(js.list[[1]])[1],js.list[[1]]$pred.m,1)]
  js.list[[2]] <- js.list[[2]][,new.aug := rtruncnorm(dim(js.list[[2]])[1],-Inf,0,js.list[[2]]$pred.m,1)]
  js.list[[3]] <- js.list[[3]][,new.aug := rtruncnorm(dim(js.list[[3]])[1],0,Inf,js.list[[3]]$pred.m,1)]
  js.list <- rbindlist(js.list)
  setkeyv(js.list,c("id","vars"))
  js <- reshape2::acast(js.list,id ~ vars,value.var = 'new.aug')
  return(as.matrix(js))
}
yy <- augment.vote.data(omega = omega, lambda = lambda, alpha = alpha, ystar = ystar, n = n, p = p)
#Generate values for bottom part of the model
#Start by getting distances between docs to get a rough clustering
#Set number of clusters to p/50
#If there are any docs that form their own cluster, reduce number of clusters
dist.docs <- dist(t(mstar), method = "euclidian", diag = T)
#Pass dist to hierarchical clustering routine
doc.groups <- hclust(dist.docs, method = "ward.D")
doc.groups <- cutree(doc.groups, k = round(p/25))
#Generate aa to have structure according to doc.groups
aa <- matrix(ncol = L, nrow = p)
for(i in 1:length(unique(doc.groups))){
  doc.base <- matrix(ncol = L, rep(rnorm(L,0,1),sum(doc.groups == sort(unique(doc.groups))[i])), byrow = T) 
  doc.vals <- matrix(ncol = L, rnorm(L*sum(doc.groups == sort(unique(doc.groups))[i]),0,.01))
  doc.base <- doc.base + doc.vals
  aa[doc.groups == sort(unique(doc.groups))[i]] <- doc.base
}
bb <- matrix(ncol = L, nrow = d, rnorm(d*L,0,1))
bb.norm <- matrix(ncol = L, rep(sqrt(1 + apply(bb^2,1,sum)),L))
bb <- bb/bb.norm
rr <- matrix(ncol = L, nrow = p, rbinom(L*p,1,.2))
aa <- aa*rr
sig <- matrix(ncol = L, nrow = d, runif(d*L,.1,1))
eta.l <- as.matrix(rbeta(L,1,1))
#Given A and B, write utility function to augment m
#Normalize over row
#THIS HAS TO BE FASTER!!!!!!
augment.text.data <- function(aa,bb,mstar,p,d){
  pred.mat <- bb%*%t(aa)
  aug.mat <- function(j){
    pm <- pred.mat[j,]
    ms <- mstar[j,]
    pmms <- data.table(cbind(seq(1,p),pm,ms))
    setnames(pmms,c("id","pm","ms"))
    setorder(pmms,"ms")
    new.augs <- c()
    um <- sort(unique(pmms$ms))
    min.vals <- -Inf
    for(i in 1:length(um)){
      c.pmms <- pmms[ms == um[i]]
      nas <- rtruncnorm(dim(c.pmms)[1], a = min.vals, b = Inf, mean = c.pmms$pm, sd = 1)
      new.augs <- c(new.augs,nas)
      min.vals <- max(new.augs)
    }
    pmms <- pmms[,new.augs := new.augs]
    setorder(pmms,"id")
    new.augs <- pmms$new.augs
    new.augs <- (new.augs - mean(new.augs))/sd(new.augs)
    return(new.augs)
  }
  aug.mat <- Vectorize(aug.mat,"j")
  new.mm <- sapply(1:d,aug.mat)
  new.mm <- t(new.mm)
  #cl <- makeCluster(4)
  #registerDoParallel(cl)
  #new.mm <- foreach(j = 1:p, .combine = "rbind",.packages = c("data.table","truncnorm")) %dopar% aug.mat(j)
  #stopCluster(cl)
  return(new.mm)
}
mm <- augment.text.data(aa = aa, bb = bb, mstar = mstar, p = p, d = d)
#Calculate the distance matrix for the initial A
kappa <- make.norm.dist.mat(aa)
#Finally, calculate pi from pi.star and kappa
pis <- kappa%*%pi.star
###Given the set of starting values, start building mcmc functions
samp.z.lam <- function(lambda,zz,alpha,omega,pis,n,p,K){
  gam <- c()
  for(i in 1:K){
    gam[i] <- (t(as.matrix(omega[,i]))%*%as.matrix(omega[,i])) + gamma.k[i]
  }
  new.lambda <- lambda
  new.zz <- zz
  for(i in 1:p){
    for(j in 1:K){
      ljk <- as.matrix(lambda[i,])
      ljk[j,1] <- 0
      eh <- t(as.matrix(yy[i,])) - (t(ljk) %*% t(omega)) + t(as.matrix(rep(alpha[i,1],n)))
      mu <- (1/gam[j])*(t(as.matrix(omega[,j]))%*%t(eh))
      pi.sum <- sum(pis[,j])
      log.rat <- .5*(log(gamma.k[j]) - log(gam[j]) + (gam[j]*(mu^2)) + log(pi.sum) - log(p - pi.sum + 1))
      if(log.rat > 10){
        zjk <- 1
      }else{
        if(log.rat < -10){
          zjk <- 0
        }else{
          z.prob <- exp(log.rat)/(1 + exp(log.rat))
          zjk <- rbinom(1,1,z.prob)
        }
      }
      if(zjk == 1){
        lamjk <- rnorm(1,mu,sqrt(1/gam[j]))
      }else{
        lamjk <- 0
      }
      new.lambda[i,j] <- lamjk
      new.zz[i,j] <- zjk
    }
    #print(i)
  }
  out.list <- list()
  out.list[[1]] <- new.lambda
  out.list[[2]] <- new.zz
  return(out.list)
}
lam.out <- samp.z.lam(lambda = lambda, zz = zz, alpha = alpha, omega = omega, pis = pis, n = n, p = p, K = K)
lambda <- lam.out[[1]]
zz <- lam.out[[2]]
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
omega <- sample.omega(yy = yy, lambda = lambda, alpha = alpha, n = n, K = K)
sample.gamma.k <- function(lambda,K){
  mk <- apply(sign(abs(lambda)),2,sum)
  sls <- apply(lambda^2,2,sum)
  new.gamma.k <- c()
  for(i in 1:K){
    new.gamma.k[i] <- rgamma(1,mk[i],sls[i])
  }
  return(new.gamma.k)
}
gamma.k <- sample.gamma.k(lambda = lambda, K = K)
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
alpha <- sample.alpha(yy = yy, lambda = lambda, omega = omega, p = p)
#Let n.near be the number of nearest neighbors to look over in this function
sample.pi.star <- function(pi,pi.star,kappa,lambda,p,K,eta.k,n.near){
  #Start by finding the n.near nearest neighbors in kappa for each obs
  near.neigh <- matrix(ncol = n.near, nrow = p)
  near.kappa <- matrix(ncol = n.near, nrow = p)
  new.pi.star <- pi.star
  new.pis <- pis
  for(i in 1:p){
    high.vals <- head(sort(kappa[i,], decreasing = T),n.near + 1)
    high.vals <- high.vals[-1]
    for(j in 1:n.near){
      near.neigh[i,j] <- which(kappa[i,] == high.vals[j])
    }
    near.kappa[i,] <- high.vals
  }
  zz <- sign(abs(lambda))
  local.z <- array(dim = c(p,K,n.near))
  for(i in 1:p){
    for(j in 1:K){
      local.z[i,j,] <- zz[near.neigh[i,],j]
    }
  }
  sum.local.z <- apply(local.z,c(1,2),sum)
  for(i in 1:p){
    #for(j in 1:K){
    #Start by pulling 20 draws from the proposal distribution
    ft <- eta.k + sum.local.z[i,]
    st <- (1 - eta.k) + (n.near - sum.local.z[i])
    props <- matrix(ncol = K, nrow = 1000)
    for(j in 1:K){
      props[,j] <- rbeta(1000,ft[j],st[j])
    }
    props[props > .999] <- .999
    props[props < .001] <- .001
    #props <- c(props,.1,.01,.3)
    #Generate current values of pi based on pi.star
    new.pis <- kappa%*%new.pi.star
    #lz is vector of nearest neighor z vals for dim k
    #slz is the sum of local z for current i and j
    #ps is the current pi.star
    #prop is the value of the proposal
    #fps is the vector of jth pi values for the nearest neighbors to i
    #ks is the vector of kappas for the nearest neighbors
    eval.mh.rat <- function(j){
      lz <- local.z[i,j,]
      slz <- sum.local.z[i,j]
      ps <- pi.star[i,j]
      prop <- props[,j] 
      fps <- new.pis[near.neigh[i,],j] 
      ks <- near.kappa[i,]
      fp <- slz*(log(ps) - log(prop))
      sp <- (n.near - slz)*(log(1 - ps) - log(1 - prop))
      tp1 <- lz*log(1 + ((ks*(prop - ps))/fps))
      tp2 <- (1 - lz)*log(1 - ((ks*(prop - ps))/(1 - fps)))
      tp <- sum(tp1 + tp2)
      return(fp + sp + tp)
    }
    prop.mh.rats <- sapply(1:K,eval.mh.rat)
    winners <- apply(prop.mh.rats,2,which.max)
    win.mat <- matrix(ncol = 3, nrow = K)
    for(j in 1:K){
      win.mat[j,1] <- new.pi.star[i,j]
      win.mat[j,2] <- props[winners[j],j]
      win.mat[j,3] <- prop.mh.rats[winners[j],j]
    }
    new.samps <- c()
    for(j in 1:K){
      if(win.mat[j,3] >= 0){
        new.samps[j] <- win.mat[j,2]
      }else{
        pp <- exp(win.mat[j,3])
        ss <- rbinom(1,1,pp)
        if(ss == 0){
          new.samps[j] <- win.mat[j,1]
        }
        if(ss == 1){
          new.samps[j] <- win.mat[j,2]
        }
      }
    }
    new.pi.star[i,] <- new.samps
    #print(i)
  }
  for(i in 1:p){
    for(j in 1:K){
      if(new.pi.star[i,j] > .999){
        new.pi.star[i,j] <- .999
      }
      if(new.pi.star[i,j] < .001){
        new.pi.star[i,j] <- .001
      }
    }
  }
  return(new.pi.star)
}
pi.star <- sample.pi.star(pi = pis, pi.star = pi.star, kappa = kappa, lambda = lambda, p = p, K = K, eta.k = eta.k, n.near = p/10)
pis <- kappa%*%pi.star
sample.eta.k <- function(eta.k,pi.star,p,K){
  mm <- pi.star/(1 - pi.star)
  mm <- log(mm)
  mm <- apply(mm,2,sum)
  mm <- -mm
  new.eta.k <- c()
  for(i in 1:K){
    uk <- runif(1,0,eta.k[i])
    vk <- runif(1,0,(1 - eta.k[i]))
    wk <- runif(1,0,(sin(pi*eta.k[i])^p))
    invwk <- asin((wk)^(1/p))/pi
    min.nk <- max(uk,min(invwk,1-invwk))
    max.nk <- min(1 - vk, max(invwk,1-invwk))
    q1 <- as.brob(-mm[i]*min.nk)
    q2 <- as.brob(-mm[i]*max.nk)
    q <- exp(q1) - exp(q2)
    nx <- runif(1,0,1)
    new.draw <- (-log(exp(q1) - (as.brob(q*nx))))/(mm[i])
    new.eta.k[i] <- new.draw
  }
  return(as.matrix(new.eta.k))
}
eta.k <- sample.eta.k(eta.k = eta.k, pi.star = pi.star, p = p, K = K)
#Set theta to reasonable value, for starts set theta to .8
theta <- .8
approx.ll <- function(new.aa,zz,pi.star){
  new.dists <- make.norm.dist.mat(new.aa)
  new.pis <- new.dists%*%pi.star
  new.probs <- (zz*new.pis) + ((1 - zz)*(1 - new.pis))
  return(sum(log(new.probs)))
}
sample.aa <- function(aa,bb,rr,mm,zz,pi.star,p,L){
  #try an element by element slice sampler
  #set.seed(1234)
  new.aa <- aa
  for(i in 1:L){
    mdj <- mm - (bb[,-i]%*%t(rr[,-i]*new.aa[,-i]))
    b2 <- apply(bb^2,2,sum)[i]
    uc.var <- 1/(1 + (rr[,i]*b2))
    bl <- matrix(ncol = p,rep(bb[,i],p))
    mbl <- bl*mdj
    mbl <- apply(mbl,2,sum)
    uc.mean <- uc.var*mbl*rr[,i]
    curr.vals <- new.aa[,i]
    curr.ll <- dnorm(curr.vals, uc.mean, sqrt(uc.var), log = T)
    aux.uc <- curr.ll - rexp(p,1)
    range.p <- uc.mean + sqrt(-uc.var*((2*aux.uc) + log(2*pi*uc.var)))
    range.m <- uc.mean - sqrt(-uc.var*((2*aux.uc) + log(2*pi*uc.var)))
    range.l <- c()
    range.u <- c()
    for(j in 1:p){
      range.l[j] <- min(range.p[j],range.m[j])
      range.u[j] <- max(range.p[j],range.m[j])
    }
    curr.bn <- approx.ll(new.aa = new.aa, zz = zz, pi.star = pi.star)
    aux.bn <- curr.bn - rexp(1,1)
    st <- 0
    cc <- 0
    while(st == 0){
      cc <- cc + 1
      naa <- new.aa
      props <- c()
      for(j in 1:p){
        props[j] <- runif(1,range.l[j],range.u[j])
      }
      naa[,i] <- props
      eval.naa <- approx.ll(new.aa = naa, zz = zz, pi.star = pi.star)
      if(eval.naa > aux.bn){
        new.aa[,i] <- props
        st <- 1
      }else{
        if(cc == 100){
          st <- 1
        }
        diffs <- sign(props - new.aa[,i])
        for(j in 1:p){
          if(diffs[j] == -1){
            range.l[j] <- props[j]
          }else{
            range.u[j] <- props[j]
          }
        }
      }
    }
    #print(i)
  }
  return(new.aa)
}
aa <- sample.aa(aa = aa, bb = bb, rr = rr, mm = mm, zz = zz, pi.star = pi.star, p = p, L = L)
sample.bb <- function(aa, bb, rr, mm, sig, d, L){
  new.bb <- bb
  aarr <- (aa^2)*rr
  aarr.l <- apply(aarr,2,sum)
  aarr.l <- matrix(nrow = d, rep(aarr.l,d))
  sig.dl <- 1/(sig + aarr.l)
  for(i in 1:L){
    mdj <- mm - (new.bb[,-i]%*%t(rr[,-i]*aa[,-i]))
    aarr.i <- aa[,i]*rr[,i]
    bb.draw <- c()
    for(j in 1:d){
      mu.dl <- sig.dl[j,i]*sum(aarr.i*mdj[j,])
      bb.draw[j] <- rnorm(1,mu.dl,sqrt(sig.dl[j,i]))
    }
    new.bb[,i] <- bb.draw
  }
  nbb.norm <- matrix(ncol = L, rep(sqrt(1 + apply(new.bb^2,1,sum)),L))
  new.bb <- new.bb/nbb.norm
  return(new.bb)
}
bb <- sample.bb(aa = aa, bb = bb, rr = rr, mm = mm, sig = sig, d = d, L = L)
sample.sig <- function(bb,sig,d,L){
  new.sig <- sig
  bb2 <- bb^2
  for(i in 1:d){
    for(j in 1:L){
      new.sig[i,j] <- rgamma(1,.5,.5*bb2[i,j])
    }
  }
  return(new.sig)
}
sig <- sample.sig(bb = bb, sig = sig, d = d, L = L)
sample.rr <- function(aa,bb,rr,mm,eta.l,p,L){
  new.rr <- rr
  for(j in 1:L){
    mdl <- mm - (bb[,-j]%*%t(new.rr[,-j]*aa[,-j]))
    bpb <- t(as.matrix(bb[,j]))%*%as.matrix(bb[,j])
    lel <- log(eta.l[j])
    l1mel <- log(1 - eta.l[j])
    for(i in 1:p){
      p1 <- lel - (.5*((bpb*aa[i,j]) - (2*t(as.matrix(bb[,j]))%*%as.matrix(mdl[,i])*aa[i,j])))
      p0 <- l1mel
      pr <- exp(p1)/(exp(p0) + exp(p1))
      pr <- max(pr,.001)
      pr <- min(pr,.999)
      new.rr[i,j] <- rbinom(1,1,pr)
    }
  }
  return(new.rr)
}
rr <- sample.rr(aa = aa, bb = bb, rr = rr, mm = mm, eta.l = eta.l, p = p, L = L)
sample.eta.l <- function(rr,p,L){
  new.eta.l <- eta.l
  s.rr <- apply(rr,2,sum)
  s.1r <- p - s.rr
  for(i in 1:L){
    b1 <- (1/1000) + s.rr[i]
    b2 <- (999/1000) + s.1r[i]
    new.eta.l[i] <- rbeta(1,b1,b2)
  }
  return(new.eta.l)
}
eta.l <- sample.eta.l(rr = rr, p = p, L = L)
pis <- make.norm.dist.mat(aa)%*%pi.star
###############################################################################################
for(m in 1:1000){
  yy <- augment.vote.data(omega = omega, lambda = lambda, alpha = alpha, ystar = ystar, n = n, p = p)
  mm <- augment.text.data(aa = aa, bb = bb, mstar = mstar, p = p, d = d)
  new.lz <- samp.z.lam(lambda = lambda, zz = zz, alpha = alpha, omega = omega, pis = pis, n = n, p = p, K = K)
  lambda <- new.lz[[1]]
  zz <- new.lz[[2]]
  omega <- sample.omega(yy = yy, lambda = lambda, alpha = alpha, n = n, K = K)
  gamma.k <- sample.gamma.k(lambda = lambda, K = K)
  alpha <- sample.alpha(yy = yy, lambda = lambda, omega = omega, p = p)
  pi.star <- sample.pi.star(pi = pis, pi.star = pi.star, kappa = kappa, lambda = lambda, p = p, K = K, eta.k = eta.k, n.near = p/10)
  eta.k <- sample.eta.k(eta.k = eta.k, pi.star = pi.star, p = p, K = K)
  aa <- sample.aa(aa = aa, bb = bb, rr = rr, mm = mm, zz = zz, pi.star = pi.star, p = p, L = L)
  kappa <- make.norm.dist.mat(aa)
  pis <- kappa%*%pi.star
  bb <- sample.bb(aa = aa, bb = bb, rr = rr, mm = mm, sig = sig, d = d, L = L)
  sig <- sample.sig(bb = bb, sig = sig, d = d, L = L)
  rr <- sample.rr(aa = aa, bb = bb, rr = rr, mm = mm, eta.l = eta.l, p = p , L = L)
  eta.l <- sample.eta.l(rr = rr, p = p, L = L)
  #Remove inactive features
  #Right now, remove based on eta
  rm.eta.l <- as.numeric(eta.l < .3)
  rm.eta.k <- apply(pis,2,sum)/p
  if(sum(rm.eta.l) != 0){
    rmel <- which(eta.l < .5)
    aa <- aa[,-rmel]
    bb <- bb[,-rmel]
    rr <- rr[,-rmel]
    sig <- sig[,-rmel]
    eta.l <- eta.l[-rmel]
    L <- length(eta.l)
  }
  if(sum(rm.eta.k) != 0){
    rmek <- which(eta.k < .5)
    lambda <- lambda[,-rmek]
    zz <- zz[,-rmek]
    omega <- omega[,-rmek]
    gamma.k <- gamma.k[-rmek]
    pi.star <- pi.star[,-rmek]
    eta.k <- eta.k[-rmek]
    kappa <- make.norm.dist.mat(aa)
    pis <- kappa%*%pi.star
    K <- length(eta.k)
  }
  print(list(length(eta.k),c(eta.k),length(eta.l),c(eta.l)))
}
