# demos

globalVariables(c("x", "olympic", "sd","marks","iris"))

demo1 = function(n=100, st=TRUE){ # approximately Gaussian quantile
  if(n > 2000 && st) stop("too large n")
  Y = ogi.value(ord.example(n))
  graphics::plot((1:n)/(n+1), Y, type="b")
  xs = (1:1000)/1001
  graphics::points(xs, stats::qnorm(xs), type="l", col="red")
}

demo2 = function(mar=FALSE){ # example of two-way table
  CT = matrix(c(
    2,1,1,0,0,
    8,3,3,0,0,
    0,2,1,1,1,
    0,0,0,1,1,
    0,0,0,0,1), 5, 5, byrow=TRUE)
  X = matrix(0, 0, 2)
  for(i in 1:5){
    for(j in 1:5){
      if(CT[i,j]>0){
        X = rbind(X, matrix(c(6-i,6-j), CT[i,j], 2, byrow=TRUE))
      }
    }
  }
  X0 = X
  X = as.data.frame(X0)
  X[,1] = factor(X0[,1], ordered=TRUE)
  X[,2] = factor(X0[,2], ordered=TRUE)
  ogiX = ogi(X, mar=mar)
  #	plot(c(0,6),c(0,6),type="n")
  #	text(X0[,1], X0[,2], round(ogiX$value,2))
  graphics::par(pty="s", cex=1.7, mar=c(4.5,3,1,1))
  graphics::plot(ogiX$scaled, xlim=c(-3,3), ylim=c(-3,3), xlab="Geometry", ylab="Probability")
  for(t in 1:nrow(ogiX$scaled)){
    xy = ogiX$scaled[t,]
    g = rep(sum(xy)/2, 2)
    graphics::segments(xy[1], xy[2], g[1], g[2], lty=2)
  }
  graphics::arrows(-3, -3, 3, 3)
  graphics::text(2.5, 2, "OGI/2")
  ogiX
}

demo3 = function(p=10, tau=0.5, rho=0.5){ # large weight
  # condition: -1/(p-2) < tau < 1, |rho| < sqrt((1+(p-2)tau)/(p-1))
  S = matrix(0, p, p)
  S[1,1] = 1
  for(j in 2:p) S[1,j] = S[j,1] = -rho
  for(i in 2:p){
    for(j in 2:p){
      if(i == j) S[i,j] = 1
      else S[i,j] = tau
    }
  }
  list(B=cov2biu(S), weight=cov2weight(S))
}

demo4 = function(se=FALSE){
  Xori = utils::read.csv("baseball2014.csv")
  X = Xori[,-c(1,2)]
  rownames(X) = X[,1]
  X[,c(12,22,23)] = X[,c(12,22,23)] * (-1)
  ogi(X, force=TRUE, se=se)
}

demo5 = function(q = 10){
  p = stats::pbinom(0:(q-1), q, 0.5) # binomial
  S = matrix(0, q, q)
  for(i in 1:q){
    for(j in 1:q){
      k = min(i,j)
      S[i,j] = p[k] - p[i] * p[j]
    }
  }
  methods::show(S)
  methods::show(eigen(S))
  d = cov2weight(S)
  e = sqrt(diag(S))
  graphics::plot(c(0, cumsum(e)), (c(0, cumsum(d)) - sum(d)/2)/sqrt(q), ylim=c(-3,3))
  d
}

demo6 = function(q=2){
  d2 = function(z, q=2) sqrt( (2+(q-1)*z^2-z*sqrt(4*q+(q-1)^2*z^2)) / 2 / (1-z^2) )
  d1 = function(z, q=2) sqrt(d2(z,q)^2+q-1)/q
  graphics::curve(d2(x,q) / d1(x,q), -1, 1)
}

demo7 = function(){
  X = datasets::iris[,1:4]
  graphics::plot(rowSums(my.scale(X)), ogi.value(X), col=iris[,5])
}

demo8 = function(eps=FALSE){ # decathlon 1988
  set.seed(20150620)
  # library(ade4)
  if (!requireNamespace("ade4", quietly = TRUE)) {
    stop("The package ade4 is required. Please install it.", call. = FALSE)
  }
  utils::data(olympic)
  Xori = olympic$tab
  # 100, long, shotput, high, 400, 110H, disc, pole, javelin, 1500
  names(Xori) = c("100m", "Long Jump", "Shot", "High Jump", "400m", "110mH", "Discus", "Pole Vault", "Javelin", "1500m")

  # scoring
  # http://www.decathlon2000.com/upload/file/pdf/scoringtables.pdf
  P = matrix(0, nrow(Xori), ncol(Xori))
  abct = matrix(c(
    25.4347, 18.00, 1.81, 1,  # 100m
    0.14354, 220.00, 1.40, 2, # Long Jump
    51.39, 1.50, 1.05, 3,     # Shot
    0.8465, 75.00, 1.42, 2,   # High Jump
    1.53775, 82.00, 1.81, 1,   # 400m
    5.74352, 28.50, 1.92, 1,  # 110mH
    12.91, 4.00, 1.10, 3,     # Discus
    0.2797, 100.00, 1.35, 2,  # Pole Vault
    10.14, 7.00, 1.08, 3,     # Javelin
    0.03768, 480.00, 1.85, 1  # 1500m
  ), 10, 4, byrow=TRUE)
  for(j in 1:10){
    a = abct[j,1]; b = abct[j,2]; c = abct[j,3]; t = abct[j,4]
    if(t == 1) P[,j] = a * (b - Xori[,j])^c  # track events
    if(t == 2) P[,j] = a * (100*Xori[,j] - b)^c  # jumps
    if(t == 3) P[,j] = a * (Xori[,j] - b)^c  # throws
  }
  P = floor(P)
  rowSums(P)  # score
  olympic$score  # score (a little different?)

  #	resP = ogi(my.scale(P), se=TRUE)
  resP = ogi(P, se=TRUE)
  graphics::plot(rowSums(P), resP$value, xlab="score", ylab="OGI")
  if(eps) grDevices::dev.copy2eps(file="demo8c.eps")

  X = Xori
  X[,c(1,5,6,10)] = -X[,c(1,5,6,10)]
  resX = ogi(my.scale(X), se=TRUE)

  graphics::par(cex=1.5, mar=c(5,4,2,2))
  nameX = ordered(names(X), names(X))
  graphics::plot(nameX, resX$weight, ylim=c(0, max(resX$weight + resX$w.se)), las=3, cex.axis=0.8, ylab="relative weight")
  graphics::segments(1:10,resX$weight-resX$w.se, 1:10,resX$weight+resX$w.se)
  #	arrows(1:10,resX$weight-resX$w.se, 1:10,resX$weight+resX$w.se, angle=90, length=0.03)
  #	arrows(1:10,resX$weight+resX$w.se, 1:10,resX$weight-resX$w.se, angle=90, length=0.03)
  if(eps) grDevices::dev.copy2eps(file="demo8a.eps")

  # weight of OGI
  graphics::plot(nameX, resP$weight, ylim=c(0, max(resP$weight + resP$w.se)), las=3, cex.axis=0.8, ylab="weight")
  graphics::segments(1:10,resP$weight-resP$w.se, 1:10,resP$weight+resP$w.se)
  if(eps) grDevices::dev.copy2eps(file="demo8b.eps")

  # weight of scale sum
  n = nrow(P)
  w = 1 / (apply(P, 2, sd) * sqrt((n-1) / n))
  se.loop = 1000
  wB = matrix(0, length(w), se.loop)
  n = nrow(P)
  for(i in 1:se.loop){
    PB = P[sample(1:n, n, replace=TRUE),]
    wB[,i] = 1 / (apply(PB, 2, sd) * sqrt((n-1) / n))
  }
  w.se = apply(wB, 1, sd)
  graphics::plot(nameX, w, ylim=c(0, max(w + w.se)), las=3, cex.axis=0.8, ylab="weight")
  graphics::segments(1:10, w - w.se, 1:10, w + w.se)
  if(eps) grDevices::dev.copy2eps(file="demo8d.eps")

  #	OGI = result$value; v.se = result$v.se
  #	ymax = max(OGI + v.se); ymin = min(OGI - v.se)
  #	plot(olympic$score, OGI, ylim=c(ymin,ymax), xlab="official record", ylab="OGI")
  #	segments(olympic$score, OGI-v.se, olympic$score, OGI+v.se)

  resP
}

demo8a = function(eps=FALSE){ # decathlon 1991 -- 2007
  Xori = utils::read.csv("decathlon.csv")
  #	Xori = Xori[Xori[,1] != 1991, -1]
  Xori = Xori[, -1]
  names(Xori) = c("100m", "Long Jump", "Shot", "High Jump", "400m", "110mH", "Discus", "Pole Vault", "Javelin", "1500m")

  Xori = Xori[Xori[,6] < 25, ]  # remove an outlier

  # scoring
  # http://www.decathlon2000.com/upload/file/pdf/scoringtables.pdf
  P = matrix(0, nrow(Xori), ncol(Xori))
  abct = matrix(c(
    25.4347, 18.00, 1.81, 1,  # 100m
    0.14354, 220.00, 1.40, 2, # Long Jump
    51.39, 1.50, 1.05, 3,     # Shot
    0.8465, 75.00, 1.42, 2,   # High Jump
    1.53775, 82.00, 1.81, 1,   # 400m
    5.74352, 28.50, 1.92, 1,  # 110mH
    12.91, 4.00, 1.10, 3,     # Discus
    0.2797, 100.00, 1.35, 2,  # Pole Vault
    10.14, 7.00, 1.08, 3,     # Javelin
    0.03768, 480.00, 1.85, 1  # 1500m
  ), 10, 4, byrow=TRUE)
  for(j in 1:10){
    a = abct[j,1]; b = abct[j,2]; c = abct[j,3]; t = abct[j,4]
    if(t == 1) P[,j] = a * (b - Xori[,j])^c  # track events
    if(t == 2) P[,j] = a * (100*Xori[,j] - b)^c  # jumps
    if(t == 3) P[,j] = a * (Xori[,j] - b)^c  # throws
  }
  P = floor(P)
  rowSums(P)  # score

  resP = ogi(P, se=TRUE)
  graphics::plot(rowSums(P), resP$value, xlab="score", ylab="OGI")
  if(eps) grDevices::dev.copy2eps(file="demo8c.eps")

  X = Xori
  X[,c(1,5,6,10)] = -X[,c(1,5,6,10)]
  resX = ogi(my.scale(X), se=TRUE)

  graphics::par(cex=1.5, mar=c(5,4,2,2))
  nameX = ordered(names(X), names(X))
  graphics::plot(nameX, resX$weight, ylim=c(0, max(resX$weight + resX$w.se)), las=3, cex.axis=0.8, ylab="relative weight")
  graphics::segments(1:10,resX$weight-resX$w.se, 1:10,resX$weight+resX$w.se)
  #	arrows(1:10,resX$weight-resX$w.se, 1:10,resX$weight+resX$w.se, angle=90, length=0.03)
  #	arrows(1:10,resX$weight+resX$w.se, 1:10,resX$weight-resX$w.se, angle=90, length=0.03)
  if(eps) grDevices::dev.copy2eps(file="demo8a.eps")

  # weight of OGI
  graphics::plot(nameX, resP$weight, ylim=c(0, max(resP$weight + resP$w.se)), las=3, cex.axis=0.8, ylab="weight")
  graphics::segments(1:10,resP$weight-resP$w.se, 1:10,resP$weight+resP$w.se)
  if(eps) grDevices::dev.copy2eps(file="demo8b.eps")

  # weight of scale sum
  n = nrow(P)
  w = 1 / (apply(P, 2, sd) * sqrt((n-1) / n))
  se.loop = 1000
  wB = matrix(0, length(w), se.loop)
  n = nrow(P)
  for(i in 1:se.loop){
    PB = P[sample(1:n, n, replace=TRUE),]
    wB[,i] = 1 / (apply(PB, 2, sd) * sqrt((n-1) / n))
  }
  w.se = apply(wB, 1, sd)
  graphics::plot(nameX, w, ylim=c(0, max(w + w.se)), las=3, cex.axis=0.8, ylab="weight")
  graphics::segments(1:10, w - w.se, 1:10, w + w.se)
  if(eps) grDevices::dev.copy2eps(file="demo8d.eps")

  #	OGI = result$value; v.se = result$v.se
  #	ymax = max(OGI + v.se); ymin = min(OGI - v.se)
  #	plot(olympic$score, OGI, ylim=c(ymin,ymax), xlab="official record", ylab="OGI")
  #	segments(olympic$score, OGI-v.se, olympic$score, OGI+v.se)

  resP
}

demo9 = function(){ # marks data
  # library(bnlearn)
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("The package bnlearn is required. Please install it.", call. = FALSE)
  }
  utils::data(marks)
  X = marks
  ogi(X, se=TRUE)
}

demo10 = function(){ # figure skating
  # Sochi, Ladies, Free Skating (FS)
  # http://www.isuresults.com/results/owg2014/SEG004.HTM
  FS = utils::read.table("FS.txt", sep="\t", header=TRUE)
  X = FS[,c("TES", "SS", "TR", "PE", "CH", "IN", "Ded.")]
  X[,7] = X[,7] * (-1)
  nu = c(.80, rep(.16, 5), .01)^2
  ogiX = ogi(X, nu=nu, se=TRUE)
  TSS = FS[,"TSS"]
  #	show(TSS - as.matrix(X) %*% c(1,rep(1.6,5),1))  # check
  graphics::plot(TSS, ogiX$value, type="n", ylab="OGI")
  graphics::text(TSS, ogiX$value, 1:nrow(X))
  ogiX
}

demo11 = function(nu=1/2, n=100, add=FALSE, lty=1, lwd=1){ # geostatistics
  if(nu <= 0) stop("not appropriate nu")
  rho = 1
  xi = seq(0, 1, len=n)
  Dxi = abs(outer(xi, xi, "-"))
  S = ifelse(Dxi != 0, 1/gamma(nu)/2^(nu-1) * (Dxi/rho)^nu * besselK(Dxi/rho, nu), 1) # Matern class
  w = cov2weight(S)
  if(!add) graphics::plot(xi, w, type="l", ylim=c(0, max(w)), lty=lty, xlab=expression(xi), lwd=lwd)
  if(add) graphics::points(xi, w, type="l", lty=lty, lwd=lwd)
  invisible(w)
}
if(FALSE){
  par(cex=1.4, mar=c(5,4,2,2))
  demo11(.1, lwd=2)
  graphics::text(0.7, 0.20, expression(paste(nu, "=0.1")))
  demo11(.5, add=TRUE, lty=2, lwd=2)
  graphics::text(0.9, 0.14, expression(paste(nu, "=0.5")))
  demo11(1, add=TRUE, lty=4, lwd=2)
  graphics::text(0.9, 0.10, expression(paste(nu, "=1.0")))
  dev.copy2eps(file="demo11.eps")
}

demo12 = function(){ # Times higher education
  # https://www.timeshighereducation.co.uk/world-university-rankings/2014-15/world-ranking
  Xori = utils::read.csv("univ.csv")
  w = !is.na(Xori[,7])  # missing values
  X = as.matrix(Xori[,4:9])/100
  # overall, teaching, international, industry, research, citations
  nu = c(100, 30, 7.5, 2.5, 30, 30)/100
  #	plot(X[,1], X[,2:6] %*% nu.sqrt[2:6])
  #	show(X[,1] - X[,2:6] %*% nu.sqrt[2:6])
  #	X = qnorm(as.matrix(Xori[w,4:9])/100)
  #	X = qnorm(as.matrix(Xori[,4:9])/100)
  ogiX = ogi(X[w,2:6], nu=nu[2:6], se=TRUE)
  graphics::plot(X[w,1], ogiX$value, type="n", ylab="OGI")
  graphics::text(X[w,1], ogiX$value, 1:nrow(X[w,]))
  ogiX
}

demo13 = function(){ # parallel coordinate plot
  X = datasets::USJudgeRatings
  lam = eigen(stats::cor(X))$vec[,1]
  graphics::par(mfrow=c(2,2))
  graphics::matplot(t(X), type="l", ylab="raw", col=1, lty=1)
  graphics::matplot(t(my.scale(X)), type="l", ylab="standardized", col=1, lty=1)
  graphics::matplot(t(my.scale(X) %*% diag(lam)), type="l", ylab="textile", col=1, lty=1)
  graphics::matplot(t(ogi.scale(X)), type="l", ylab="OGI", col=1, lty=1)
}

demo14 = function(){ # USJudgeRatings
  set.seed(20150620)
  X = datasets::USJudgeRatings
  p = ncol(X)
  resX = ogi(X, se=TRUE)
  graphics::par(cex=1.5, mar=c(5,4,2,2))
  nameX = ordered(names(X), names(X))
  graphics::plot(nameX, resX$weight, ylim=c(0, max(resX$weight + resX$w.se)), las=3, cex.axis=0.8, ylab="relative weight")
  graphics::segments(1:p, resX$weight-resX$w.se, 1:p, resX$weight+resX$w.se)
  #	dev.copy2eps(file="demo14.eps")
}

if(FALSE){  # garbage
  king = read.table("king_eng.txt", header=TRUE)
  king[,c(1, 5, 6, 10)] = king[,c(1, 5, 6, 10)] * (-1)

  test1 = matrix(-7/12, 3, 3)
  test1[2:3,2:3] = 0
  diag(test1) = rep(1,3)
}
