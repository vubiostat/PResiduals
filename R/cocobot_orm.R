###### cocobot.orm #####
##### correlation of probility scale residuals, but fit cumulative probability models instead of linear models
library(SparseM)
library(rms)

corTS = function(xresid, yresid,
                 xz.dl.dtheta, yz.dl.dtheta,
                 xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                 dxresid.dthetax, dyresid.dthetay,fisher=FALSE){
  
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
  
  if (fisher==TRUE){
    ####Fisher's transformation
    TS_f <- log( (1+TS)/(1-TS) )
    var.TS_f <- varTS*(2/(1-TS^2))^2
    pvalTS <- 2 * pnorm( -abs(TS_f)/sqrt(var.TS_f))
  }
  
  list(TS=TS,var.TS=varTS, pval.TS=pvalTS)
}

getCI <- function(ts, var, fisher, ci=0.95){
  if(!fisher){
    lower <- ts - abs(qnorm(0.5*(1-ci)))*sqrt(var)
    upper <- ts + abs(qnorm(0.5*(1-ci)))*sqrt(var)  
  } else {
    ts_f <- log( (1+ts)/(1-ts) )
    var_f <- var*(2/(1-ts^2))^2
    lower_f <- ts_f - abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    upper_f <- ts_f + abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    lower <- (exp(lower_f)-1)/(1+exp(lower_f))
    upper <- (exp(upper_f)-1)/(1+exp(upper_f))
  }
  return(c(lower, upper))
}


######
orm.scores <- function (y, X, link=c("logistic", "probit", "cauchit", "loglog", "cloglog")) {
  y <- as.numeric(factor(y))
  if(!is.matrix(X)) X = matrix(X, ncol=1)
  N = length(y)
  ny = length(table(y))
  na = ny - 1
  nb = ncol(X)
  npar = na + nb
  Z = outer(y, 1:ny, "<=")
  
  mod <- orm(y~X, family=link[1], x=TRUE, y=TRUE)
  
  alpha = mod$coeff[1:na]
  beta = mod$coeff[-(1:na)]
  
  dl.dtheta = matrix(0, N, npar)
  #   ## Information matrices are stored as sums over all individuals.
  #   d2l.dalpha.dalpha = matrix(0,na,na)
  #   d2l.dalpha.dbeta = matrix(0,na,nb)
  #   d2l.dbeta.dbeta = matrix(0,nb,nb)
  #   d2l.dbeta.dalpha = matrix(0,nb,na)
  
  #p0 = matrix(,N,ny)
  #dp0.dtheta = array(0,c(N,ny,npar))
  ## Cumulative probabilities
  Gamma = matrix(0,N,na)
  #dgamma.dtheta = array(0,c(N,na,npar))
  dlow.dtheta = dhi.dtheta =dpij.dtheta= matrix(, npar, N)
  
  
  for (i in 1:N) {
    z = Z[i,]  ## z has length ny
    x = X[i,]
    
    ## gamma and phi are defined as in McCullagh (1980)
    #####gamma = 1 - 1/(1 + exp(alpha + sum(beta*x))) ## gamma has length na
    ### I am using orm() notation, orm model P(Y>=j+1), gamma=P(Y<=j)=1- F(alpha+sum(beta*x))
    ### It does not matter much for probit, logit, cauchy (with symetric pdf), but it does matter for loglog and cloglog
    ### Therefore,I modify the cocobot ordinal.scores code.
    gamma <-1- mod$trans$cumprob(alpha + sum(beta*x))
    diffgamma = diff(c(gamma,1))
    invgamma = 1/gamma
    invgamma2 = invgamma^2
    invdiffgamma = 1/diffgamma
    invdiffgamma2 = invdiffgamma^2
    phi = log(gamma / diffgamma) ## phi has length na
    Gamma[i,] = gamma
    
    
    #### Some intermediate derivatives
    ## g(phi) = log(1+exp(phi))
    dg.dphi = 1 - 1/(1 + exp(phi))
    ## l is the log likelihood (6.3) in McCullagh (1980)
    dl.dphi = z[-ny] - z[-1] * dg.dphi
    t.dl.dphi = t(dl.dphi)
    
    ## dphi.dgamma is a na*na matrix with rows indexed by phi
    ## and columns indexed by gamma
    dphi.dgamma = matrix(0,na,na)
    diag(dphi.dgamma) = invgamma + invdiffgamma
    if(na > 1)
      dphi.dgamma[cbind(1:(na-1), 2:na)] = -invdiffgamma[-na]
    
    
    dgamma.base <- - mod$trans$deriv(x=alpha + sum(beta*x), f=mod$trans$cumprob(alpha + sum(beta*x)))
    
    dgamma.dalpha = diagn(dgamma.base)
    dgamma.dbeta = dgamma.base %o% x
    #dgamma.dtheta[i,,] = cbind(dgamma.dalpha, dgamma.dbeta)
    #### First derivatives of log-likelihood (score functions)    
    dl.dalpha = dl.dphi %*% dphi.dgamma %*% dgamma.dalpha
    dl.dbeta = dl.dphi %*% dphi.dgamma %*% dgamma.dbeta
    dl.dtheta[i,] = c(dl.dalpha, dl.dbeta)
    
    if (y[i] == 1) {
      dlow.dtheta[,i] <- 0
    } else {
      dlow.dtheta[,i] <- c(dgamma.dalpha[y[i]-1,], dgamma.dbeta[y[i]-1, ])
    }
    
    if (y[i] == ny) {
      dhi.dtheta[,i] <- 0
    } else {
      
      dhi.dtheta[,i] <- -c(dgamma.dalpha[y[i],], dgamma.dbeta[y[i],])
    }
    
    
    
    
    
    
    
    
    
    ######################### use orm information matrix directly
    ####### I have checked they give same results with the following calculation
    #     d2gamma.base= -mod$trans$deriv2(x=alpha + sum(beta*x), f=mod$trans$cumprob(alpha + sum(beta*x)),
    #                                     deriv=mod$trans$deriv(x=alpha + sum(beta*x), f=mod$trans$cumprob(alpha + sum(beta*x))))
    #     
    #  
    #     ##
    #     d2l.dphi.dphi = diagn(-z[-1] * dg.dphi * (1-dg.dphi))
    #     d2l.dphi.dalpha = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dalpha
    #     d2l.dphi.dbeta = d2l.dphi.dphi %*% dphi.dgamma %*% dgamma.dbeta
    #     
    #     ##
    #     d2phi.dgamma.dalpha = array(0,c(na,na,na))
    #     d2phi.dgamma.dalpha[cbind(1:na,1:na,1:na)] = (-invgamma2 + invdiffgamma2) * dgamma.base
    #     if(na > 1) {
    #       d2phi.dgamma.dalpha[cbind(1:(na-1),1:(na-1),2:na)] = -invdiffgamma2[-na] * dgamma.base[-1]
    #       d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,1:(na-1))] = -invdiffgamma2[-na] * dgamma.base[-na]
    #       d2phi.dgamma.dalpha[cbind(1:(na-1),2:na,2:na)] = invdiffgamma2[-na] * dgamma.base[-1]
    #     }
    #     
    #     ##
    #     d2phi.dgamma.dbeta = array(0,c(na,na,nb))
    #     rowdiff = matrix(0,na,na)
    #     diag(rowdiff) = 1
    #     if(na > 1) rowdiff[cbind(1:(na-1),2:na)] = -1
    #     d2phi.dgamma.dbeta.comp1 = diagn(-invdiffgamma2) %*% rowdiff %*% dgamma.dbeta
    #     d2phi.dgamma.dbeta.comp2 = diagn(-invgamma2) %*% dgamma.dbeta - d2phi.dgamma.dbeta.comp1
    #     for(j in 1:na) {
    #       d2phi.dgamma.dbeta[j,j,] = d2phi.dgamma.dbeta.comp2[j,]
    #       if(j < na)
    #         d2phi.dgamma.dbeta[j,j+1,] = d2phi.dgamma.dbeta.comp1[j,]
    #     }
    #     
    #     ##
    #     d2gamma.dalpha.dbeta = array(0,c(na,na,nb))
    #     for(j in 1:na)
    #       d2gamma.dalpha.dbeta[j,j,] = d2gamma.base[j] %o% x
    #     
    #     ##
    #     d2gamma.dbeta.dbeta = d2gamma.base %o% x %o% x
    
    
    
    
    
    #### Second derivatives of log-likelihood
    #### Since first derivative is a sum of terms each being a*b*c,
    #### second derivative is a sum of terms each being (a'*b*c+a*b'*c+a*b*c').
    
    #### d2l.dalpha.dalpha
    ## Obtain aprime.b.c
    ## Transpose first so that matrix multiplication is meaningful.
    ## Then transpose so that column is indexed by second alpha.
    #     aprime.b.c = t(crossprod(d2l.dphi.dalpha, dphi.dgamma %*% dgamma.dalpha))
    #     
    #     ## Obtain a.bprime.c
    #     ## run through the index of second alpha
    #     a.bprime.c = matrix(,na,na)
    #     for(j in 1:na)
    #       a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dalpha[,,j] %*% dgamma.dalpha
    #     
    #     ## Obtain a.b.cprime
    #     ## cprime = d2gamma.dalpha.dalpha = 0 if indices of the two alphas differ.
    #     d2gamma.dalpha.dalpha = diagn(d2gamma.base)
    #     a.b.cprime = diagn(as.vector(dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dalpha))
    #     
    #     ## summing over individuals
    #     d2l.dalpha.dalpha = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dalpha
    #     
    #     
    #     #### d2l.dalpha.dbeta
    #     aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dalpha))
    #     a.bprime.c = a.b.cprime = matrix(,na,nb)
    #     for(j in 1:nb) {
    #       a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dalpha
    #       a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dalpha.dbeta[,,j]
    #     }
    #     d2l.dalpha.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dalpha.dbeta
    #     
    #     
    # 
    #     
    #     #### d2l.dbeta.dbeta
    #     aprime.b.c = t(crossprod(d2l.dphi.dbeta, dphi.dgamma %*% dgamma.dbeta))
    #     a.bprime.c = a.b.cprime = matrix(,nb,nb)
    #     for(j in 1:nb) {
    #       a.bprime.c[,j] = t.dl.dphi %*% d2phi.dgamma.dbeta[,,j] %*% dgamma.dbeta
    #       a.b.cprime[,j] = t.dl.dphi %*% dphi.dgamma %*% d2gamma.dbeta.dbeta[,,j]
    #     }
    #     d2l.dbeta.dbeta = aprime.b.c + a.bprime.c + a.b.cprime + d2l.dbeta.dbeta 
    #     
    
    #### Derivatives of predicted probabilities
    #     p0[i,] = diff(c(0, gamma, 1))
    #     
    #     rowdiff = matrix(0,ny,na)
    #     diag(rowdiff) = 1
    #     rowdiff[cbind(2:ny,1:na)] = -1
    #     dp0.dalpha = rowdiff %*% dgamma.dalpha
    #     dp0.dbeta = rowdiff %*% dgamma.dbeta
    #     
    #     dp0.dtheta[i,,] = cbind(dp0.dalpha, dp0.dbeta)
    #     
    #     dpij.dtheta[,i] <- dp0.dtheta[i, y[i],]
    #     
    #     
    #     if (y[i] == 1) {
    #       dlow.dtheta[,i] <- 0
    #     } else {
    #       dlow.dtheta[,i] <- dgamma.dtheta[i,y[i]-1,]
    #     }
    #     
    #     if (y[i] == ny) {
    #       dhi.dtheta[,i] <- 0
    #     } else {
    #       dhi.dtheta[,i] <- -dgamma.dtheta[i,y[i],]
    #     }
    
    
  }
  
  #   d2l.dtheta.dtheta = rbind(
  #     cbind(d2l.dalpha.dalpha, d2l.dalpha.dbeta),
  #     cbind(t(d2l.dalpha.dbeta), d2l.dbeta.dbeta))
  #   
  
  low = cbind(0, Gamma)[cbind(1:N, y)]
  hi = cbind(1-Gamma, 0)[cbind(1:N, y)]
  presid <- low - hi
  dpresid.dtheta <- dlow.dtheta - dhi.dtheta
  
  result <- list(dl.dtheta=dl.dtheta,
                 d2l.dtheta.dtheta = -mod$info.matrix,
                 presid=presid,
                 dpresid.dtheta = dpresid.dtheta )
  
  #   if(latent){
  #     rij <- cbind(Gamma, 1)[cbind(1:N, y)]
  #     rij_1 <- cbind(0,Gamma)[cbind(1:N, y)]
  #     pij <- rij-rij_1
  #     
  #     #G.inverse <- switch(link[1], logistic = qlogis, probit = qnorm,
  #     #                    cloglog = qgumbel, cauchit = qcauchy)
  #     G.inverse <- mod$trans$inverse
  #     latent.resid <- rep(NA, N)
  #     
  #     inverse_fail <- FALSE 
  #     for (i in 1:N){
  #       tmp <- try(integrate(G.inverse, rij_1[i], rij[i])$value/pij[i],silent=TRUE)
  #       if (inherits(tmp,'try-error')){
  #         if (link[1] != 'cauchit')
  #           warning("Cannot compute latent variable residual.")
  #         else
  #           warning("Cannot compute latent variable residual with link function cauchit.")
  #         inverse_fail <- TRUE
  #         break
  #       } else {
  #         latent.resid[i] <- tmp
  #       }
  #     }
  #     
  #     dlatent.dtheta = matrix(, npar, N)
  #     
  #     if (!inverse_fail){
  #       ### To compute dlatent.dtheta (need dgamma.dtheta and dp0.dtheta from ordinal scores)
  #       
  #       
  #       drij_1.dtheta <- dlow.dtheta
  #       drij.dtheta <- - dhi.dtheta
  #       for(i in 1:N) {  
  #         
  #         
  #         if (y[i] == 1) {
  #           dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
  #             G.inverse(rij[i])*drij.dtheta[,i] - 0 )
  #         } else if(y[i] == ny){
  #           dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
  #             0 - G.inverse(rij_1[i])*drij_1.dtheta[,i] )
  #         } else
  #           dlatent.dtheta[,i] <- -latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
  #             G.inverse(rij[i])*drij.dtheta[,i] - G.inverse(rij_1[i])*drij_1.dtheta[,i])
  #       }
  #     }
  #     
  #     result <- list(dl.dtheta=dl.dtheta,
  #                    d2l.dtheta.dtheta = -mod$info.matrix,
  #                    presid=presid,
  #                    dpresid.dtheta = dpresid.dtheta,
  #                    latent.resid=latent.resid,
  #                    dlatent.dtheta=dlatent.dtheta
  #     )
  #     
  #     
  #   }
  
  
  return(result)
  
}


cocobot.orm <- function(formula, data, link.x=c("logistic", "probit", "cauchit", "loglog", "cloglog"),
                        link.y=c("logistic", "probit", "cauchit", "loglog", "cloglog"),
                        subset, na.action=getOption('na.action'), 
                        fisher=FALSE,conf.int=0.95) {
  
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
  
  ## return value
  ans <- list(
    TS=list(),
    fisher=fisher,
    conf.int=conf.int,
    data.points=data.points
  )
  
  score.xz <- orm.scores(y=model.response(mx), X=mxz, link=link.x)
  score.yz <- orm.scores(y=model.response(my), X=myz, link=link.y)
  ts = corTS(score.xz$presid, score.yz$presid,
             score.xz$dl.dtheta, score.yz$dl.dtheta,
             as.matrix(score.xz$d2l.dtheta.dtheta), as.matrix(score.yz$d2l.dtheta.dtheta),
             score.xz$dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
  ts.label = "PResid vs. PResid"
  
  ans$TS$TS <-  list( ts=ts$TS, var=ts$var.TS, pval=ts$pval.TS,
                      label = ts.label)
  ans <- structure(ans, class="cocobot")
  
  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }
  
  ans
  
  
}

# ##### example
# ##### same result with cobot, if both X and Y are ordinal variables
# generate.data = function(alphax, betax, alphay, betay, eta, N) {
#   z = rnorm(N,0,1)
#   x = y = numeric(N)
#   
#   ## px is an N x length(alphax) matrix.
#   ## Each row has the TRUE cummulative probabilities for each subject.
#   px = (1 + exp(- outer(alphax, betax*z, "+"))) ^ (-1)
#   aa = runif(N)
#   for(i in 1:N)
#     x[i] = sum(aa[i] > px[,i])
#   x = as.numeric(as.factor(x))
#   ## x = x+1 may have category gaps if there are small probability categories.
#   
#   py = (1 + exp(- outer(alphay, betay*z+eta[x], "+"))) ^ (-1)
#   aa = runif(N)
#   for(i in 1:N)
#     y[i] = sum(aa[i] > py[,i])
#   y = as.numeric(as.factor(y))
#   ## y = y+1 may have category gaps if there are small probability categories.
#   
#   return(list(x=x, y=y, z=z))
# }


# N = 500
# alphay = c(-1, 0, 1)
# betay = -.5
# alphax = c(-1, 0, 1, 2)
# betax = 1
# eta = rep(0,5)
# data = generate.data(alphax, betax, alphay, betay, eta, N)
# cobot(as.factor(x)|as.factor(y) ~z, data=data)
# cocobot.orm(x|y~z, data=data)
# cobot(as.factor(x)|as.factor(y) ~z, data=data, link="probit")
# cocobot.orm(x|y~z, data=data, link.x="probit", link.y="probit") 
# 
# 
# ####### when X and Y are continuous
# set.seed(1)
# n <- 500
# z <- rnorm(n, 0, 1)
# x.latent <- rnorm(n, z+0.5, 1)
# x <-  exp(x.latent) + x.latent^3
# y.latent <- rnorm(n, x.latent+0.5*z)
# y <- exp(y.latent)
# #cobot(as.factor(x) |as.factor(y) ~z)
# Sys.time()
# cocobot.orm(x|y~z, link.x="logistic", link.y="logistic")
# Sys.time()
# cocobot.orm(x|y~z, link.x="probit", link.y="probit")
