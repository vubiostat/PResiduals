poisson.scores <- function(y, X, ...){
  N <-  length(y)
  mod <- glm(y~X, family=poisson(), ...)
  #dl.dtheta <- estfun(mod)
  dl.dtheta <- cbind(y- mod$fitted.value, (y- mod$fitted.value)*X)
  #d2l.dtheta.dtheta <- -solve(bread(mod))*N
  d2l.dtheta.dtheta <- - crossprod(cbind(1, X)*sqrt(mod$fitted))
  
  presid <- ppois(y-1, mod$fitted.values) + ppois(y, mod$fitted.values) -1
  
  #dpresid.dlambda <- 2*(ppois(y-2, mod$fitted.values) - ppois(y-1, mod$fitted.values)) + dpois(y-1, mod$fitted.values) - dpois(y, mod$fitted.values)
  #### can be further simplified as
  dpresid.dlambda <- -  (dpois(y-1, mod$fitted.values) + dpois(y, mod$fitted.values))
  dlambda.dtheta <- mod$fitted.values * cbind(1, X)
  dpresid.dtheta <- t(dpresid.dlambda * dlambda.dtheta)
  
  pearson.resid <- (y - mod$fitted.values)/sqrt(mod$fitted.values)
  dpearson.resid.dlambda <- - 0.5 * ( mod$fitted.values^(-1/2) + mod$fitted.values^(-3/2)*y)
  dpearson.resid.dtheta <- t(dpearson.resid.dlambda * dlambda.dtheta)
  
  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       presid = presid,
       dpresid.dtheta = dpresid.dtheta,
       pearson.resid = pearson.resid,
       dpearson.resid.dtheta = dpearson.resid.dtheta)
}



nb.scores <- function(y, X, ...){
  N <-  length(y)
  #mod <- glm.nb(y~X)
  mod <- glm.nb(y~X,  ...)
  #### score for beta's
  dl.dbeta <- (y - mod$fitted.values)/(1 + mod$fitted.values/mod$theta) * cbind(1, X)
  #### score for 1/size  (size = theta in glm.nb)
  dl.dsize <- unlist(lapply(y, FUN=function(y) ifelse(y==0, 0, sum(1/(mod$theta + 0:(y-1)))))) - log( 1 + mod$fitted.values/mod$theta) + (mod$fitted.values - y)/(mod$theta + mod$fitted.values)
  dl.dtheta <- cbind(dl.dbeta, dl.dsize)
  
  #### fisher information for beta's
  d2l.dbeta.dbeta <- crossprod(cbind(1, X)*sqrt(mod$fitted.values/(1+mod$fitted.values/mod$theta)))
  #### fisher information for theta
  d2l.dsize.dsize <- sum( unlist(lapply(y, FUN=function(y) ifelse(y==0, 0, sum(1/(mod$theta + 0:(y-1))^2)))) - mod$fitted.values/mod$theta/(mod$fitted.values + mod$theta))
  d2l.dtheta.dtheta <- rbind(cbind(d2l.dbeta.dbeta,0),0)
  d2l.dtheta.dtheta[dim(d2l.dtheta.dtheta)[1] , dim(d2l.dtheta.dtheta)[2] ] <- d2l.dsize.dsize
  presid <- pnbinom(y-1, mu=mod$fitted.values, size=mod$theta ) + pnbinom(y, mu=mod$fitted.values, size=mod$theta) -1
  pearson.resid <- (y - mod$fitted.values)/sqrt(mod$fitted.values + mod$fitted.values^2/mod$theta)
  dmu.dbeta <- mod$fitted.values * cbind(1, X)
  dpresid.dmu <- dnbinom(y-1, size=mod$theta, mu=mod$fitted.values)/(mod$fitted.value + mod$theta) - (dnbinom(y-1, size=mod$theta, mu=mod$fitted.values) + dnbinom(y, size=mod$theta, mu=mod$fitted.values))*(mod$theta + y)/(mod$theta + mod$fitted.values)
  dpresid.dbeta <- t(dpresid.dmu * dmu.dbeta)
  df.dsize <- function(k, mu, size) {
    dnbinom(k, size=size, mu=mu)*((mu-k)/(mu+size) +
                                     log(size) -log(mu+size) +
                                     unlist(lapply(k, FUN=function(k) ifelse(k==0, 0, sum(1/(size + (0:(k-1))))))) )                                  
  }
  
  ### dpresid.dsize = 2 dF(y-1).dsize + df(y).dsize = 2 sum(df(y-1).dsize) + df(y).dsize
  #tmp <-  cbind(y, mod$fitted.values)
  #tmp2 <- apply(tmp, 1, FUN=function(x) ifelse(x[1]==0, 0, 2*sum(df.dsize(k=0:(x[1]-1), mu=x[2], size=mod$theta))))
 
  dpresid.dsize <- df.dsize(y, mod$fitted, mod$theta)  +
    apply(cbind(y, mod$fitted.values), 1, FUN=function(x) ifelse(x[1]==0, 0, 2*sum(df.dsize(k=0:(x[1]-1), mu=x[2], size=mod$theta))))
  dpresid.dtheta <- rbind(dpresid.dbeta, dpresid.dsize)
  
  dpearson.resid.dmu <- -0.5 * (y * (mod$fitted.values + mod$fitted.values^2/mod$theta)^(-3/2)*(1+2*mod$fitted.values/mod$theta) +
                                  (1/mod$fitted.values + 1/mod$theta)^(-3/2)/mod$fitted.values^2)
  dpearson.resid.dbeta <- t(dpearson.resid.dmu * dmu.dbeta)
  dpearson.resid.dsize <- 0.5*(y-mod$fitted.values)*(mod$fitted.values + mod$fitted.values^2/mod$theta)^(-3/2)*mod$fitted.values^2/mod$theta^2
  dpearson.resid.dtheta <- rbind(dpearson.resid.dbeta, dpearson.resid.dsize)
  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       presid = presid,
       dpresid.dtheta = dpresid.dtheta,
       pearson.resid = pearson.resid,
       dpearson.resid.dtheta = dpearson.resid.dtheta)
  

}

##### same function as in cocobot
#corTS = function(xresid, yresid,
#                 xz.dl.dtheta, yz.dl.dtheta,
#                 xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
#                 dxresid.dthetax, dyresid.dthetay,fisher=FALSE){
#  
#  TS = cor(xresid, yresid)
#  
#  xresid2 = xresid^2
#  yresid2 = yresid^2
#  xbyyresid = xresid * yresid
#  mean.xresid = mean(xresid)
#  mean.yresid = mean(yresid)
#  mean.xbyyresid = mean(xbyyresid)
#  
#  bigphi = cbind(xz.dl.dtheta,
#                 yz.dl.dtheta,
#                 mean.xresid - xresid,
#                 mean.yresid - yresid,
#                 mean.xbyyresid - xbyyresid,
#                 mean(xresid2)-xresid2,
#                 mean(yresid2)-yresid2,
#                 0)
#  
#  npar.xz = dim(xz.dl.dtheta)[2]
#  npar.yz = dim(yz.dl.dtheta)[2]
#  Ntheta = npar.xz + npar.yz + 6
#  N = dim(xz.dl.dtheta)[1]
#  
#  A = matrix(0,Ntheta,Ntheta)
#  A[1:npar.xz, 1:npar.xz] = xz.d2l.dtheta.dtheta
#  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = yz.d2l.dtheta.dtheta
#  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)
#  
#  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
#                     c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
#                     c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
#                     c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
#                     c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
#  
#  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial
#  
#  ## TS also equals numTS / sqrt(varprod) = numTS * revsvp
#  numTS = mean.xbyyresid - mean.xresid * mean.yresid
#  var.xresid = mean(xresid2) - mean.xresid^2
#  var.yresid = mean(yresid2) - mean.yresid^2
#  varprod = var.xresid * var.yresid
#  revsvp = 1/sqrt(varprod)
#  dTS.dvarprod = numTS * (-0.5) * revsvp^3
#  
#  smallpartial = N *
#    c(-mean.yresid * revsvp + dTS.dvarprod * (-2*mean.xresid*var.yresid),
#      -mean.xresid * revsvp + dTS.dvarprod * (-2*mean.yresid*var.xresid),
#      revsvp,
#      dTS.dvarprod * var.yresid,
#      dTS.dvarprod * var.xresid)
#  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
#  
#  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
#  
#  SS = solve(A, t(bigphi))
#  var.theta = tcrossprod (SS, SS)
#  varTS = var.theta[Ntheta, Ntheta]
#  pvalTS = 2 * pnorm( -abs(TS)/sqrt(varTS))
#  
#  if (fisher==TRUE){
#    ####Fisher's transformation
#    TS_f <- log( (1+TS)/(1-TS) )
#    var.TS_f <- varTS*(2/(1-TS^2))^2
#    pvalTS <- 2 * pnorm( -abs(TS_f)/sqrt(var.TS_f))
#  }
#  
#  list(TS=TS,var.TS=varTS, pval.TS=pvalTS)
#}

### formula x|y ~z, x is ordinal and y is continous
## family argument is added to specifies the family (poisson or negative binomial) for count data (new for countbot)
## emp argument is deleted compared with cocobot
countbot <- function(formula, data, link=c("logit", "probit", "cloglog", "cauchit"),
                     family=c("poisson", "negative binomial"),
                     subset, na.action=getOption('na.action'), 
                     fisher=FALSE, conf.int=0.95) {
  
  
  # Construct the model frames for x ~ z and y ~ z
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  # We set xlev to a benign non-value in the call so that it won't get partially matched
  # to any variable in the formula. For instance a variable named 'x' could possibly get
  # bound to xlev, which is not what we want.
  mf$xlev <- integer(0) 
  mf[[1L]] <- as.name("model.frame")
  
  
  mx <- my <- mf
  
  # NOTE: we add the opposite variable to each model frame call so that
  # subsetting occurs correctly. Later we strip them off.
  mx[["formula"]] <- Fx
  yName <- all.vars(Fy[[2]])[1]
  mx[[yName]] <- Fy[[2]]
  
  my[["formula"]] <- Fy
  xName <- all.vars(Fx[[2]])[1]
  my[[xName]] <- Fx[[2]]
  
  mx <- eval(mx, parent.frame())
  mx[[paste('(',yName,')',sep='')]] <- NULL
  
  my <- eval(my, parent.frame())
  my[[paste('(',xName,')',sep='')]] <- NULL
  
  data.points <- nrow(mx)
  
  if (!is.factor(mx[[1]])){
    warning("Coercing ",names(mx)[1]," to factor. Check the ordering of categories.")
    mx[[1]] <- as.factor(mx[[1]])
  }
  
  if (is.factor(my[[1]])){
    stop(names(my)[1]," cannot be a factor.")
  }
  
  # Construct the model matrix z
  mxz <- model.matrix(attr(mx,'terms'),mx) 
  zzint <- match("(Intercept)", colnames(mxz), nomatch = 0L)
  if(zzint > 0L) {
    mxz <- mxz[, -zzint, drop = FALSE]
  }
  
  myz <- model.matrix(attr(my,'terms'),my) 
  zzint <- match("(Intercept)", colnames(myz), nomatch = 0L)
  if(zzint > 0L) {
    myz <- myz[, -zzint, drop = FALSE]
  }
  
  score.xz <- ordinal.scores(mx, mxz,method=link)
  if (family[1]=="poisson")
    score.yz <- poisson.scores(y=model.response(my), X=myz)
  else if (family[1]=="negative binomial")
    score.yz <- nb.scores(y=model.response(my), X=myz)
  else stop("family has to be 'poisson' or 'negative binomial'")
  
  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]
  xx = as.integer(model.response(mx))
  
  nx = length(table(xx))
  
  N = length(xx)
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]
  
  xz.presid <- low.x - hi.x
  xz.dpresid.dtheta <- score.xz$dlow.dtheta - score.xz$dhi.dtheta
  
  ## return value
  ans <- list(
    TS=list(),
    fisher=fisher,
    conf.int=conf.int,
    data.points=data.points
  )
  
  tb = corTS(xz.presid, score.yz$presid,
             score.xz$dl.dtheta, score.yz$dl.dtheta,
             score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
             xz.dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
  tb.label = "PResid vs. PResid (assume normality)" 
  
  ans$TS$TB <- 
    list( 
      ts=tb$TS, var=tb$var.TS, pval=tb$pval.TS,
      label = tb.label
    )
  
  rij <- cbind(score.xz$Gamma, 1)[cbind(1:N, xx)]
  rij_1 <- cbind(0,score.xz$Gamma)[cbind(1:N, xx)]
  pij <- rij-rij_1
  
  G.inverse <- switch(link[1], logit = qlogis, probit = qnorm,
                      cloglog = qgumbel, cauchit = qcauchy)
  xz.latent.resid <- rep(NA, N)
  
  inverse_fail <- FALSE 
  for (i in 1:N){
    tmp <- try(integrate(G.inverse, rij_1[i], rij[i])$value/pij[i],silent=TRUE)
    if (inherits(tmp,'try-error')){
      if (link[1] != 'cauchit')
        warning("Cannot compute latent variable residual.")
      else
        warning("Cannot compute latent variable residual with link function cauchit.")
      inverse_fail <- TRUE
      break
    } else {
      xz.latent.resid[i] <- tmp
    }
  }
  
  if (!inverse_fail){
    ### To compute dlatent.dtheta (need dgamma.dtheta and dp0.dtheta from ordinal scores)
    xz.dlatent.dtheta = dpij.dtheta = matrix(, npar.xz, N)
    
    drij_1.dtheta <- score.xz$dlow.dtheta
    drij.dtheta <- -score.xz$dhi.dtheta
    for(i in 1:N) {  
      dpij.dtheta[,i] <- score.xz$dp0.dtheta[i, xx[i],]
      
      if (xx[i] == 1) {
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - 0 )
      } else if(xx[i] == nx){
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          0 - G.inverse(rij_1[i])*drij_1.dtheta[,i] )
      } else
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - G.inverse(rij_1[i])*drij_1.dtheta[,i])
    }
    
    
    ### latent.resid vs pearson resid
    tc <- corTS(xz.latent.resid, score.yz$pearson.resid,
                score.xz$dl.dtheta, score.yz$dl.dtheta,
                score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                xz.dlatent.dtheta, score.yz$dpearson.resid.dtheta, fisher)
    
    ans$TS$TC <- 
      list( 
        ts=tc$TS, var=tc$var.TS, pval=tc$pval.TS,
        label = 'Latent.resid vs. Pearson.resid'
      )
  }
  ans <- structure(ans, class="cocobot")
  
  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }
  
  ans
  
}


#### example
## generate count by ordinal data
generate.data3 = function(alphax, betax, alphay, betay, eta, N) {
  z = rnorm(N,0,1)
  x = y = numeric(N)
  
  ## px is an N x length(alphax) matrix.
  ## Each row has the TRUE cummulative probabilities for  each subject.
  px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
  aa = runif(N)
  for(i in 1:N)
    x[i] = sum(aa[i] > px[,i])
  x = as.numeric(as.factor(x))
  y = rpois(N, exp(outer(alphay, betay*z+eta[x], "+")))
  
  return(list(x=as.factor(x), y=y, z=z))
}


set.seed(13)
alphax = c(-1, 0, 1, 2)
betax = 1
alphay = 1
betay = -.5

#eta = rep(0, 5)
eta = c(1:5)/20 
N = 200
data <- generate.data3(alphax, betax, alphay, betay, eta, N)

#### check for cocobot
cocobot(x|y~z, data=data)
countbot(x|y~z, data=data, fisher=TRUE)
countbot(x|y~z, data=data, family="negative binomial")