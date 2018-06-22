#Augment Text Data
#Assumes complete data set - i.e. no NAs
#Doesn't really make sense to have missingness in this matrix
#Computationall intense, but data.table optimized
#Arguments - aa is matrix of lower level latent vars, bb is lower level factor loadings 
#mstar is the matrix of observed counts for each token: number of tokens x number of docs
#p is number of docs
#d is number of tokens
augment.text.data <- function(aa,bb,mstar,p,d){
  pred.mat <- bb%*%t(aa)
  mstar.table <- apply(mstar,1,table)
  mstar.table <- lapply(mstar.table,cumsum)
  normalize.vec <- function(x){
    x1 <- x/p
    x2 <- c(0,head(x/p,length(x) - 1))
    x1 <- qnorm(x1)
    x2 <- qnorm(x2)
    xx <- cbind(as.numeric(names(x1)),x2,x1)
    colnames(xx) <- NULL
    rownames(xx) <- NULL
    return(xx)
  }
  mstar.table <- lapply(mstar.table,normalize.vec)
  for(i in 1:d){
    mstar.table[[i]] <- data.table(cbind(c(i),mstar.table[[i]]))
  }
  mstar.table <- rbindlist(mstar.table)
  setnames(mstar.table,c("row.id","obs","lower.b","upper.b"))
  mstar.melt <- data.table(reshape2::melt(mstar))
  setnames(mstar.melt,c("row.id","col.id","obs"))
  pred.melt <- data.table(reshape2::melt(pred.mat))
  setnames(pred.melt,c("row.id","col.id","pred"))
  setkeyv(mstar.melt,c("row.id","col.id"))
  setkeyv(pred.melt,c("row.id","col.id"))
  big.melt <- pred.melt[mstar.melt]
  rm(mstar.melt,pred.melt)
  gc()
  #1,30 -2.545530
  setkeyv(big.melt,c("row.id","obs"))
  setkeyv(mstar.table,c("row.id","obs"))
  big.join.tab <- big.melt[mstar.table]
  rm(big.melt)
  gc()
  big.join.tab <- big.join.tab[,new.val := rtruncnorm(1,lower.b,upper.b,pred,1)]
  out.tab <- reshape2::acast(big.join.tab,row.id ~ col.id, value.var = "new.val")
  rownames(out.tab) <- NULL
  colnames(out.tab) <- NULL
  return(out.tab)
}