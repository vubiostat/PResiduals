###glm()
#' @export
presid.glm <- function(mod, ...) {
    ##gaussian.emp = (2 * rank(residuals(mod))-1-length(y))/length(y)  #need to figure out how to incorporate this in this function
    switch(mod$family$family ,
           poisson = 2 * ppois(mod$y, mod$fitted.values) - dpois(mod$y, mod$fitted.values) - 1,
           binomial = mod$y - mod$fitted.values,
           gaussian = 2 * pnorm((mod$y - mod$fitted.values)/sqrt(summary(mod)$dispersion)) - 1,
           stop("Unhandled family", mod$family$family))
}

###lm()
#' @export
presid.lm <- function(mod, ...) {
  presid = 2 * pnorm((model.response(mod$model) - mod$fitted.values)/summary(mod)$sigma) - 1
  ##presid.emp = (2 * rank(residuals(mod))-1-length(y))/length(y)
  
  presid
}

###negative binomial
#' @export
presid.negbin <- function(mod, ...) {
  pnbinom(mod$y-1, mu=mod$fitted.values, size=mod$theta ) + pnbinom(mod$y, mu=mod$fitted.values, size=mod$theta) -1
}


###polr
#' @export
presid.polr <- function(mod, ...) {
    pfun <- switch(mod$method,
                   logistic = plogis,
                   probit = pnorm, 
                   loglog = pgumbel,
                   cloglog = pGumbel,
                   cauchit = pcauchy)
    n <- length(mod$lp)
    q <- length(mod$zeta)
    cumpr <- cbind(0, matrix(pfun(matrix(mod$zeta, n, q, byrow = TRUE) - mod$lp),, q), 1)
    y <- as.integer(model.response(mod$model))
    lo <- cumpr[cbind(seq_len(n), y)]
    hi <- cumpr[cbind(seq_len(n), y+1L)] - 1
    lo - hi
}

###coxph()
#' @export
presid.coxph <- function(mod, ...) {
    time <- mod$y[,1]
    delta <- mod$y[,2]
    
    1-exp(residuals(mod)-delta)-delta*exp(residuals(mod)-delta)
}


###survreg()
#' @export
presid.survreg <- function(mod, ...){
    time <- mod$y[,1]
    delta <- mod$y[,2]
    
    switch(mod$dist,
           weibull = pweibull(time, shape=1/summary(mod)$scale, scale=exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pweibull(time, shape=1/summary(mod)$scale, scale=exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
           
           exponential = pexp(time, rate=1/exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pexp(time, rate=1/exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
           
           gaussian = pnorm(time, mean=mod$linear.predictors, sd=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(pnorm(time, mean=mod$linear.predictors, sd=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
           
           logistic = plogis(time, location=mod$linear.predictors, scale=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(plogis(time, location=mod$linear.predictors, scale=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
         
           loglogistic = pllogis(time, shape=summary(mod)$scale, scale=exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pllogis(time, shape=summary(mod)$scale, scale=exp(mod$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
         
           lognormal = plnorm(time, meanlog=mod$linear.predictors, sdlog=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(plnorm(time, meanlog=mod$linear.predictors, sdlog=summary(mod)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
           stop("Unhandled dist", mod$dist))
}

#' @export
presid.default <- function(mod, ...) {
    stop("Unhandeled model type")
}


#' Probability-scale Residual
#'
#' \code{presid} Calculates the probability-scale residual for various model types.
#' probability-scale residual is \eqn{P(Y* > y) - P(Y* < y)} where \eqn{y} is the observed
#' outcome and \eqn{Y*} is a random variable from the fitted distribution.
#'
#' @param mod The model object for which the probability-scale residual is calculated
#' @param ... Aditional arguements passed to methods
#' @return The probablity scale residual for the model
#' @references Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480
#' @references Shepherd BE, Li C, Lin Q.  Probability-scale residuals for continuous,
#' discrete, and censored data.  Submitted.
#' @author Charles Dupont \email{charles.dupont@@vanderbilt.edu}
#' @author Chun Li \email{chun.li@@vanderbilt.edu}
#' @author Bryan Shepherd \email{bryan.shepherd@@vanderbilt.edu}
#' @importFrom actuar pllogis
#' @importFrom stats plnorm pnorm pexp pweibull plogis pnbinom
#' @export
#' @examples
#' library(survival)
#' library(stats)
#' 
#' set.seed(100)
#' time <- sample(1:60, size=100, replace=TRUE)
#' delta <- sample(c(0,1), size=100, prob=c(0.2,0.8), replace=TRUE)
#' X <- round(rnorm(100, mean=36, sd=7), 0)
#'
#' mod <- survreg(Surv(time, delta) ~ X, dist="weibull")
#' summary(presid(mod))
#'
#' #Example for coxph
#' mod <- coxph(Surv(time, delta) ~ X)
#' summary(presid(mod))
#' ##################
#'
#'
#' #Example for negative binomial
#' library(MASS)
#' n <- 2000
#' beta0 <- 0
#' beta1 <- 1
#' X <- rnorm(n)
#' mu <- beta0 + beta1*X
#' phi <- 3
#' y <- rnbinom(n, mu=exp(mu), size=phi)
#'
#' mod <- glm.nb(y ~ X)
#' summary(presid(mod))
presid <- function(mod, ...) {
    UseMethod('presid', mod)
}

