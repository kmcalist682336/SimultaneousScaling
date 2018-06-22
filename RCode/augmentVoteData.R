#Augment vote data
#Assumes that vote data is 0, 1, or NA
#1 is a yes, 0 is a no, NA is no vote
#NA assumed to MAR (address this someday when I have tenure)
#Sets augment break at 0
#Arguments - omega is the upper level matrix of latent scores (n x K), lambda is upper level factor loadings (p x K)
#Alpha is the vector of item level intercepts (p vector)
#ystar is the matrix of observed votes (p times n)
#n is the number of voters
#p is the number of items/documents/votes
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