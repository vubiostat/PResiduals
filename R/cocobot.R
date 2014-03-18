ord.resid <- function(y, x){
  if (!is.matrix(x)) x = matrix(x, ncol=1)
  ny = length(table(y))
  nay = ny - 1
  N = length(y)
  assign("EE", 0, pos=1)
  suppressWarnings(tryCatch(assign("mod", lrm(y ~ x, tol=1e-50, maxit=100)),
                            error = function(x) assign("EE", 1, pos=1)))
  ## Obtain predicted cumulative probabilities.  Each subject is a column.
  ## Because lrm() fits y>=c instead of y<=c, the cumulative probabilities
  ## are 1- 1/[1+exp(-alpha-beta*z)] = 1/[1+exp(alpha+beta*z)]
  py = (1 + exp(outer(mod$coeff[1:nay],
                      as.vector(x%*%mod$coeff[-(1:nay)]), "+"))) ^ (-1)
  if(EE == 1) return(presid=NULL )
  
  low.y = rbind(0,py)[cbind(y,1:N)]
  hi.y = 1 - rbind(py,1)[cbind(y,1:N)]
  
  rij <- rbind(py,1)[cbind(y,1:N)]
  rij_1  <- rbind(0,py)[cbind(y,1:N)]
  pij = rij -rij_1
  G.inverse <- function(p) log(p / (1 - p))
  
  latent.resid <- rep(NA, N)
  for (i in 1:N){
    latent.resid[i] <- integrate(G.inverse, rij_1[i], rij[i])$value/pij[i]
  }
  
  
  list(presid = low.y - hi.y, latent.resid=latent.resid, rij=rij, rij_1=rij_1)
}


lm.score = function(y, X){
  N = length(y)  
  mod = lm(y~X)
  smod = summary(mod)
  ## bread = [1/N sum (- partial phi)]^-1 
  ##       = [- 1/N  d2l. dtheta. dtheta]^-1
  d2l.dtheta.dtheta = - solve(bread(mod))*N
  dl.dtheta = estfun(mod)
  resid = smod$residuals
  presid = 2*pnorm((y - mod$fitted.values)/smod$sigma) -1
  dresid.dtheta = t(cbind(-1, -X))
  dpresid.dtheta = t(cbind(-2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma,
                           -2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma *
                             X))
  
  f.y<-density(resid)
  fy.ry <- NULL
  presid.k <- NULL
  for (i in 1:length(resid)){
    fy.ry[i] <- f.y$y[which(abs(f.y$x-resid[i])==min(abs(f.y$x-resid[i])))]
    presid.k[i] <- sum(resid<resid[i])/length(resid) - sum(resid>resid[i])/length(resid)
  }
  dpresid.dtheta.k <- t(cbind(-2*fy.ry,
                              -2*fy.ry*X))
  list(mod = mod, 
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = resid,
       dresid.dtheta = dresid.dtheta,
       presid = presid,
       presid.k= presid.k,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = dpresid.dtheta.k)
}

###copy from cobot, add dgammaij_1.dtheta


ordinal.dresid.dtheta = function(ordinal.scores, y){

  N = length(y)
  ny = length(table(y))
  yy = y
  npar = dim(ordinal.scores$dl.dtheta)[2]
  dgamma.dtheta = ordinal.scores$dgamma.dtheta
  dlow.dtheta = dhi.dtheta = matrix(, npar, N)
  for(i in 1:N) {
    if (yy[i] == 1) {
      dlow.dtheta[,i] <- 0
    } else {
      dlow.dtheta[,i] <- dgamma.dtheta[i,yy[i]-1,]
    }
    
    
    if (yy[i] == ny) {
      dhi.dtheta[,i] <- 0
    } else {
      dhi.dtheta[,i] <- -dgamma.dtheta[i,yy[i],]
    }
  }
  dresid.dtheta = dlow.dtheta - dhi.dtheta
  
}

dlatent.dtheta <- function(ord.resid, ordinal.scores,y){
  G.inverse <- function(p) log(p / (1 - p))
  N = length(y)
  ny = length(table(y))
  yy = y
  npar = dim(ordinal.scores$dl.dtheta)[2]
  dgamma.dtheta = ordinal.scores$dgamma.dtheta
  dp0.dtheta = ordinal.scores$dp0.dtheta
  
  latent.resid <- ord.resid$latent.resid
  rij <- ord.resid$rij
  rij_1 <- ord.resid$rij_1
  pij <- rij-rij_1
  
  dlatent.dtheta = drij.dtheta = drij_1.dtheta = dpij.dtheta = matrix(, npar, N)
 
  for(i in 1:N) {
    if (yy[i] == 1) {
      drij_1.dtheta[,i] <- 0
    } else {
      drij_1.dtheta[,i] <- dgamma.dtheta[i,yy[i]-1,]
    }
       
    if (yy[i] == ny) {
      drij.dtheta[,i] <- 0
    } else {
      drij.dtheta[,i] <- dgamma.dtheta[i,yy[i],]
    }
    
    dpij.dtheta[,i] <- dp0.dtheta[i, yy[i],]
    
    if (yy[i] == 1) {
      dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
        G.inverse(rij[i])*drij.dtheta[,i] - 0 )
    } else if(yy[i] == ny){
      dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
        0 - G.inverse(rij_1[i])*drij_1.dtheta[,i] )
    } else
      dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
        G.inverse(rij[i])*drij.dtheta[,i] - G.inverse(rij_1[i])*drij_1.dtheta[,i])
  }
  
  return(dlatent.dtheta)
}

cor.TS = function(xresid, yresid,
                  xz.dl.dtheta, yz.dl.dtheta,
                  xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                  dxresid.dthetax, dyresid.dthetay){
  
  TS = cor(xresid, yresid)
  
  xresid2 = xresid^2
  yresid2 = yresid^2
  xbyyresid = xresid * yresid
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  
  bigphi = cbind(xz.dl.dtheta,
                 yz.dl.dtheta,
                 mean.xresid - xresid,
                 mean.yresid - yresid,
                 mean.xbyyresid - xbyyresid,
                 mean(xresid2)-xresid2,
                 mean(yresid2)-yresid2,
                 0)
  
  npar.xz = dim(xz.dl.dtheta)[2]
  npar.yz = dim(yz.dl.dtheta)[2]
  Ntheta = npar.xz + npar.yz + 6
  N = dim(xz.dl.dtheta)[1]
  
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = xz.d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = yz.d2l.dtheta.dtheta
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)
  
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
                     c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
                     c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial
  
  ## TS also equals numTS / sqrt(varprod) = numTS * revsvp
  numTS = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dTS.dvarprod = numTS * (-0.5) * revsvp^3
  
  smallpartial = N *
    c(-mean.yresid * revsvp + dTS.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dTS.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dTS.dvarprod * var.yresid,
      dTS.dvarprod * var.xresid)
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod (SS, SS)
  varTS = var.theta[Ntheta, Ntheta]
  pvalTS = 2 * pnorm( -abs(TS)/sqrt(varTS))
  ####Fisher's transformation
  TS_f <- log( (1+TS)/(1-TS) )
  var.TS_f <- varTS*(2/(1-TS^2))^2
  pvalTS_f <- 2 * pnorm( -abs(TS_f)/sqrt(var.TS_f))
  
  list(TS=TS,var.TS=varTS, pval.TS=pvalTS, TS.f=TS_f, var.TS_f=var.TS_f, pval.TS_f=pvalTS_f)
}

cocobot <- function(formula, link=c("logit", "probit", "cloglog", "cauchit"),
                  link.x=link,
                  link.y=link,
                  data, subset) {
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.fail
  mf[[1L]] <- as.name("model.frame")

  mx <- my <- mf

  mx[["formula"]] <- Fx
  my[["formula"]] <- Fy

  mx <- eval(mx, parent.frame())
  my <- eval(my, parent.frame())

  xz.score = ordinal.scores(mx,method=link.x)
  yz.score = lm.score(data$y, data$z)
  xz.resid = ord.resid(as.integer(data$x), data$z)
  xresid = xz.resid$presid
  xlatent.resid = xz.resid$latent.resid
  
  ypresid = yz.score$presid
  ypresid.k = yz.score$presid.k
  yresid = yz.score$resid
  
  xz.dresid.dtheta = ordinal.dresid.dtheta (xz.score, as.integer(data$x))
  xz.dlatent.resid.dtheta = dlatent.dtheta(ord.resid=xz.resid,ordinal.scores=xz.score, y=as.integer(data$x))
  yz.dpresid.dtheta = yz.score$dpresid.dtheta
  yz.dresid.dtheta = yz.score$dresid.dtheta
  
  ### presid vs obs-exp
#  ts1 <-  cor.TS(xresid, yresid,
#                 xz.score$dl.dtheta, yz.score$dl.dtheta,
#                 xz.score$d2l.dtheta.dtheta, yz.score$d2l.dtheta.dtheta,
#                 xz.dresid.dtheta, yz.dresid.dtheta)
  
  ### presid vs presid (use pdf of normal)
  ts2 = cor.TS(xresid, ypresid,
               xz.score$dl.dtheta, yz.score$dl.dtheta,
               xz.score$d2l.dtheta.dtheta, yz.score$d2l.dtheta.dtheta,
               xz.dresid.dtheta, yz.dpresid.dtheta)

  ### presid vs presid (emprical)
  ts2.k = cor.TS(xresid, ypresid.k,
                 xz.score$dl.dtheta, yz.score$dl.dtheta,
                 xz.score$d2l.dtheta.dtheta, yz.score$d2l.dtheta.dtheta,
                 xz.dresid.dtheta, yz.score$dpresid.dtheta.k)

  structure(list(T2=ts2$TS, varT2=ts2$var.TS, pvalT2=ts2$pval.TS,T2.E=ts2.k$TS, varT2.E=ts2.k$var.TS,pvalT2.E=ts2.k$pval.TS), class="cocobot")
}



