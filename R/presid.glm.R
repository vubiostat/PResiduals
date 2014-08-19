###glm()
#' @export
presid.glm <- function(object, emp=FALSE, ...) {
    ##gaussian.emp = (2 * rank(residuals(object))-1-length(y))/length(y)  #need to figure out how to incorporate this in this function
    switch(object$family$family ,
           poisson = 2 * ppois(object$y, object$fitted.values) - dpois(object$y, object$fitted.values) - 1,
           binomial = object$y - object$fitted.values,
           gaussian = if(emp) (2 * rank(residuals(object)) - 1 - length(y)) / length(y) else 2 * pnorm((object$y - object$fitted.values)/sqrt(summary(object)$dispersion)) - 1,
           stop("Unhandled family", object$family$family))
}

###lm()
#' @export
presid.lm <- function(object, emp=FALSE, ...) {
    y <- model.response(object$model)
  if(emp) {
      (2 * rank(residuals(object)) - 1 - length(y)) / length(y)
  } else {
      2 * pnorm((y - object$fitted.values)/summary(object)$sigma) - 1
  }
}

###negative binomial
#' @export
presid.negbin <- function(object, ...) {
  pnbinom(object$y-1, mu=object$fitted.values, size=object$theta ) + pnbinom(object$y, mu=object$fitted.values, size=object$theta) -1
}


###polr
#' @export
presid.polr <- function(object, ...) {
    pfun <- switch(object$method,
                   logistic = plogis,
                   probit = pnorm, 
                   loglog = pgumbel,
                   cloglog = pGumbel,
                   cauchit = pcauchy)
    n <- length(object$lp)
    q <- length(object$zeta)
    cumpr <- cbind(0, matrix(pfun(matrix(object$zeta, n, q, byrow = TRUE) - object$lp),, q), 1)
    y <- as.integer(model.response(object$model))
    lo <- cumpr[cbind(seq_len(n), y)]
    hi <- cumpr[cbind(seq_len(n), y+1L)] - 1
    lo - hi
}

###coxph()
#' @export
presid.coxph <- function(object, ...) {
    time <- object$y[,1]
    delta <- object$y[,2]
    
    1-exp(residuals(object)-delta)-delta*exp(residuals(object)-delta)
}


###survreg()
#' @export
presid.survreg <- function(object, ...){
    time <- object$y[,1]
    delta <- object$y[,2]
    
    switch(object$dist,
           weibull = pweibull(time, shape=1/summary(object)$scale, scale=exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pweibull(time, shape=1/summary(object)$scale, scale=exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
           
           exponential = pexp(time, rate=1/exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pexp(time, rate=1/exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
           
           gaussian = pnorm(time, mean=object$linear.predictors, sd=summary(object)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(pnorm(time, mean=object$linear.predictors, sd=summary(object)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
           
           logistic = plogis(time, location=object$linear.predictors, scale=summary(object)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(plogis(time, location=object$linear.predictors, scale=summary(object)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
         
           loglogistic = pllogis(time, shape=summary(object)$scale, scale=exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE)
           + delta*(pllogis(time, shape=summary(object)$scale, scale=exp(object$linear.predictors), lower.tail=TRUE, log.p=FALSE) - 1),
         
           lognormal = plnorm(time, meanlog=object$linear.predictors, sdlog=summary(object)$scale, lower.tail=TRUE, log.p=FALSE)
           + delta*(plnorm(time, meanlog=object$linear.predictors, sdlog=summary(object)$scale, lower.tail=TRUE, log.p=FALSE) - 1),
           stop("Unhandled dist", object$dist))
}

#' @export
presid.default <- function(object, ...) {
    stop("Unhandled model type")
}


#' Probability-scale Residual
#'
#' \code{presid} Calculates the probability-scale residual for various model function objects.
#' 
#' Probability-scale residual is \eqn{P(Y* < y) - P(Y* > y)} where \eqn{y} is the observed
#' outcome and \eqn{Y*} is a random variable from the fitted distribution.
#'
#' @param object The model object for which the probability-scale residual is calculated
#' @param ... Additional arguements passed to methods
#' @return The probability scale residual for the model
#' @references Shepherd BE, Li C, Liu Q.  Probability-scale residuals for continuous,
#' discrete, and censored data.  Submitted.
#' @references Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480
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
presid <- function(object, ...) {
    UseMethod('presid', object)
}

