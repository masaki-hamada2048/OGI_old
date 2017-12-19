#' Weight Vectors of the Bi-unit Canonical Form
#'
#' \code{cov2weight(S)} returns the numeric vector in which the diagonal
#' elements of the matrix \eqn{D} are arranged, where \eqn{DSD} is the bi-unit
#' canonical form of \eqn{S}.
#'
#' @param S Covariance matrix, especially it is positive semi-definite.
#' @param Dvec Numeric vector of initial values of iteration.
#' @param nu Numeric vector of subjective importance. It determines the
#'   importance of each of the variates.
#' @param tol Numeric number of tolerance. If the minimum eigenvalue of \code{S}
#'   is less than \code{tol}, \code{S} is considered not to be positive
#'   definite.
#' @param force Logical: if force=FALSE, \code{S} should be strictly positive
#'   definite. Default: FALSE.
#' @return Numeric vector of diagonal elements of \eqn{D}, which appears in the
#'   bi-unit canonical form \eqn{DSD} of \eqn{S}.
#' @examples
#' S = matrix(0, 5, 5)
#' S[1,1] = 1
#' for(j in 2:5) S[1,j] = S[j,1] = -0.5
#' for(i in 2:5){
#'   for(j in 2:5){
#'     if(i == j) S[i,j] = 1
#'     else S[i,j] = 0.5
#'   }
#' }
#' weight=cov2weight(S)
#' weight
#' @export
cov2weight = function(S, Dvec = rep(1, nrow(S)), nu = rep(1, nrow(S)), tol = 1e-06, force = FALSE){
  mineig = min(eigen(S)$val)
  if(mineig < tol){
    if(!force) stop("S is not positive definite.")
    if(mineig < -tol) stop("S is not positive semi-definite.")
    if(cov.check(S)$status == 0) stop("ker(S) has a non-negative vector.")
  }
  if(any(nu<=0)){
    stop("nu is not positive.")
  }
  nms = rownames(S)
  E = sqrt(diag(S))
  S = diag(1/E) %*% S %*% diag(1/E)
  p = nrow(S)
  ep = Inf
  while(ep > tol){
    D2 = Dvec
    for(i in 1:p){ # coordinate descent
      b = sum(S[i,-i]*D2[-i])
      r = sqrt(b^2 + 4*S[i,i]*nu[i])
      D2[i] = ifelse(b>0, 2*S[i,i]*nu[i]/(b+r), (-b+r)/2/S[i,i])
    }
    ep = max(abs(D2 - Dvec))
    Dvec = D2
  }
  Dvec = Dvec / E
  names(Dvec) = nms
  Dvec
}

#' Bi-unit Canonical Form
#'
#' \code{cov2biu(S)} returns the bi-unit canonical form of \code{S}.
#'
#' @param S Covariance matrix, especially it is positive semi-definite.
#' @param nu Numeric vector of subjective importance. It determines the
#'   importance of each of the variates.
#' @param force Logical: if force=FALSE, \code{S} should be strictly positive
#'   definite. Default: FALSE.
#' @param detail Logical: if detail=TRUE, it returns the list of the bi-unit
#'   form and the weight vectors. Default: FALSE.
#' @return Numeric matrix of the bi-unit canonical form \eqn{DSD} of \eqn{S}.
#' @examples
#' S = matrix(0, 5, 5)
#' S[1,1] = 1
#' for(j in 2:5) S[1,j] = S[j,1] = -0.5
#' for(i in 2:5){
#'   for(j in 2:5){
#'     if(i == j) S[i,j] = 1
#'     else S[i,j] = 0.5
#'   }
#' }
#' B=cov2biu(S)
#' B
#' @export
cov2biu = function(S, nu = rep(1, nrow(S)), force = FALSE, detail=FALSE){ # covariance matrix to bi-unit matrix
  w = cov2weight(S, nu=nu, force=force)
  B = diag(w) %*% S %*% diag(w)
  dimnames(B) = dimnames(S)
  if(detail) return(list(B=B, weight=w))
  B
}

#' Objective General Index
#'
#' \code{ogi(X)} returns the objective general index (OGI) of the covariance
#' matrix \code{S} of \code{X}.
#'
#' Consider a data matrix of \eqn{n} individuals with \eqn{p} variates. The
#' objective general index (OGI) is a general index that combines the \eqn{p}
#' variates into a univariate index in order to rank the \eqn{n} individuals.
#' The OGI is always positively correlated with each of the variates. For more
#' details, see the references.
#'
#' @param X Numeric or ordered matrix.
#' @param se Logical: if se=TRUE, it additionally computes \code{w.se} and
#'   \code{v.se} by bootstrap. Default: FALSE.
#' @param force Logical: if force=FALSE, \code{S} should be strictly positive
#'   definite. Default: FALSE.
#' @param se.loop Iteration number in bootstrap for computation of standard
#'   error.
#' @param nu Numeric vector of subjective importance. It determines the
#'   importance of each column of \code{X}.
#' @param center Logical: if center=TRUE, \code{ogi(X)$Z} is centered.
#'   Default:TRUE.
#' @param mar Logical: if mar=TRUE, each of ordered categorical variates of
#'   \code{X} (if exists) is marginally converted into a numeric vector in
#'   advance by the univariate OGI quantification. If mar=FALSE, the
#'   simultaneous OGI quantification is applied. Default:FALSE.
#' @return
#' \item{value}{The objective general index (OGI).}
#' \item{X}{The input matrix \code{X}.}
#' \item{scaled}{The product of \code{Z \%*\% diag(weight)}, where \code{Z} and
#'   \code{weight} are as follows.}
#' \item{Z}{Numerical matrix converted from \code{X}. If center = TRUE, it is centered.}
#' \item{weight}{The output of \code{\link{cov2weight}(S, nu=nu, force=force)},
#'   where \code{S} is the covariance matrix of \code{X}. }
#' \item{rel.weight}{The product of \code{weight * sqrt(diag(S))}, where \code{S}
#'   is the covariance matrix of \code{X}.}
#' \item{biu}{The bi-unit canonical form of the covariance matrix of \code{X}.}
#' \item{idx}{Numeric vector. If \code{X} has ordered categorical variates,
#'   \code{idx} has (number of levels) -1 number of indexes.}
#' \item{w.se}{If requested, \code{w.se} is numeric vector of the standard error
#'   of \code{weight}. It is calculated by bootstrap.}
#' \item{v.se}{If requested, \code{v.se} is numeric vector of the standard error
#'   of \code{value}. It is calculated by bootstrap.}
#' @references  Sei, T. (2016). An objective general index for multivariate
#'   ordered data, Journal of Multivariate Analysis, 147, 247-264.
#'   \url{http://www.sciencedirect.com/science/article/pii/S0047259X16000269}
#' @examples
#' CT = matrix(c(
#' 2,1,1,0,0,
#' 8,3,3,0,0,
#' 0,2,1,1,1,
#' 0,0,0,1,1,
#' 0,0,0,0,1), 5, 5, byrow=TRUE)
#' X = matrix(0, 0, 2)
#' for(i in 1:5){
#'   for(j in 1:5){
#'     if(CT[i,j]>0){
#'       X = rbind(X, matrix(c(6-i,6-j), CT[i,j], 2, byrow=TRUE))
#'     }
#'   }
#' }
#' X0 = X
#' X = as.data.frame(X0)
#' X[,1] = factor(X0[,1], ordered=TRUE)
#' X[,2] = factor(X0[,2], ordered=TRUE)
#' ogiX = ogi(X)
#' par(pty="s", cex=1.7, mar=c(4.5,3,1,1))
#' plot(ogiX$scaled, xlim=c(-3,3), ylim=c(-3,3), xlab="Geometry", ylab="Probability")
#' for(t in 1:nrow(ogiX$scaled)){
#'   xy = ogiX$scaled[t,]
#'   g = rep(sum(xy)/2, 2)
#'   segments(xy[1], xy[2], g[1], g[2], lty=2)
#' }
#' arrows(-3, -3, 3, 3)
#' text(2.5, 2, "OGI/2")
#' ogiX
#'
#'
#' f = ordered(1:10)
#' f[sample(1:10, 20, replace=TRUE)]
#' Y = ogi(f)$value
#' plot((1:10)/(10+1), Y, type="b")
#' xs = (1:1000)/1001
#' points(xs, qnorm(xs), type="l", col="red")
#'
#'
#' X = USJudgeRatings
#' ogiX = ogi(X)
#' nameX = ordered(names(X), names(X))
#' plot(nameX, ogiX$weight, las=3, cex.axis=0.8, ylim=c(0,1.2), ylab="weight")
#' @export
ogi = function(X, se=FALSE, force=FALSE, se.loop=1000, nu=rep(1, ncol(X)), center=TRUE, mar=FALSE){ # The main function
  if(!is.data.frame(X)) X = as.data.frame(X)
  p = ncol(X)
  if(mar){
    for(i in 1:p){
      if(is.ordered(X[,i]))
        r = ogi(X[,i,drop=FALSE], force=force, mar=FALSE)
      X[,i] = r$scaled[,1]
    }
  }
  Z = ogi.numeric(X)
  idx = attributes(Z)$index
  nu.ori = nu
  nu = rep(1, ncol(Z))
  for(i in 1:p){
    nu[idx==i] = nu.ori[i]/sum(idx==i)
  }
  Z = scale(Z, scale=FALSE, center=center)
  if(center){
    S = my.cov(Z)
  }else{
    S = t(Z) %*% Z / nrow(Z)
  }
  w = cov2weight(S, nu=nu, force=force)
  rw = w * sqrt(diag(S))
  Zw = Z %*% diag(w)
  scaled = matrix(0, nrow(Z), p)
  for(i in 1:p){
    scaled[,i] = rowSums(Zw[,idx==i,drop=FALSE])
  }
  OGI = rowSums(scaled)
  B = diag(w) %*% S %*% diag(w)
  result = list(value = OGI, X = X, scaled = scaled, Z = Z, weight = w, rel.weight = rw, biu = B, idx = idx)

  if(se){ # standard error by bootstrap
    n = nrow(X)
    wB = matrix(0, length(w), se.loop)
    for(i in 1:se.loop){
      ZB = Z[sample(1:n, n, replace=TRUE),]
      ogiB = ogi(ZB, se=FALSE, force=TRUE)
      wB[,i] = ogiB$weight
    }
    w.se = apply(wB, 1, stats::sd)
    v.se = sqrt(Z^2 %*% w.se^2)
    ymax = max(OGI + v.se); ymin = min(OGI - v.se)
    graphics::plot(rowSums(my.scale(Z)), OGI, ylim=c(ymin,ymax), xlab="scaled sum", ylab="OGI")
    #		text(rowSums(my.scale(Z)), OGI, 1:n)
    graphics::segments(rowSums(my.scale(Z)), OGI-v.se, rowSums(my.scale(Z)), OGI+v.se)
    result$w.se = w.se
    result$v.se = v.se
  }

  result
}

ogi.value = function(X){
  ogi(X)$value
}

ogi.scale = function(X){
  ogi(X)$scaled
}

ogi.weight = function(X){
  ogi(X)$weight
}

ogi.rel.weight = function(X){
  ogi(X)$rel.weight
}

ogi.biu = function(X){
  ogi(X)$biu
}

# ToDo
# predict.ogi = function(ogiX, Y){
#   X = ogiX$X
#   is.matrix(Y)
# }

# ToDo
# table2factor = function(tb){
#   if(!is.table(tb)) tb = as.table(tb)
#   d = length(dim(tb))
#   as.data.frame(tb)
# }
