# util

my.cov = function(X){
  n = nrow(X)
  stats::cov(X) * (n-1) / n
}

my.scale = function(X){
  n = nrow(X)
  scale(X) / sqrt((n-1) / n)
}

cov.check = function(S){
  # library(lpSolve)
  p = nrow(S)
  result = lpSolve::lp("min", rep(0,p), rbind(S, rep(1,p)), rep("=",p+1), c(rep(0,p),1))
  result
}

ord2num = function(v){ # ordered to numeric matrix
  if(!is.ordered(v)) stop("class not match.")
  n = length(v)
  q = nlevels(v)
  if(q <= 1) stop("too few levels.")
  lv = levels(v)
  Y = matrix(0, n, q-1) # contrast
  for(j in 2:q){
    w = (v == lv[j])
    Y[w,1:(j-1)] = 1
  }
  Y
}

ord.example = function(q, n=NULL){
  f = ordered(1:q)
  if(is.null(n)) return(f)
  f[sample(1:q, n, replace=TRUE)]
}

ogi.numeric = function(X){ # data.frame to matrix
  n = nrow(X)
  p = ncol(X)
  Z = matrix(NA, n, 0)
  idx = numeric(0)
  for(i in 1:p){
    if(is.numeric(X[,i])){
      Z = cbind(Z, X[,i])
      idx = c(idx, i)
      next
    }
    if(is.ordered(X[,i])){
      Y = ord2num(X[,i])
      Z = cbind(Z, Y)
      idx = c(idx, rep(i, ncol(Y)))
      next
    }
    stop("class not match.")
  }
  attr(Z, "index") = idx
  Z
}
