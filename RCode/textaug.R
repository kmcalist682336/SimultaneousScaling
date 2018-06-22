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
  new.mm <- sapply(1:p)
  new.mm <- t(new.mm)
  #cl <- makeCluster(4)
  #registerDoParallel(cl)
  #new.mm <- foreach(j = 1:p, .combine = "rbind",.packages = c("data.table","truncnorm")) %dopar% aug.mat(j)
  #stopCluster(cl)
  return(new.mm)
}