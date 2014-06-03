#' cobot class print method
#' @param x cobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
#' @method print cobot
#' @author Charles Dupont

print.cobot <- function(x, ...) {
  y <- matrix(nrow=length(x$TS),ncol=5)
  dims <- character(length(x$TS))
  for (i in 1:length(x$TS)){
    y[i,] <- c(x$TS[[i]]$ts, sqrt(x$TS[[i]]$var),x$TS[[i]]$pval,x$TS[[i]]$lower,x$TS[[i]]$upper)
    dims[i] <- x$TS[[i]]$label
  }
  dimnames(y) <- list(dims,c('est','stderr','p','lower CI','upper CI'))
  invisible(print(y,...))
  cat('Fisher Transform:',x$fisher,'\n')
  cat('Confidence Interval:',x$conf.int,'\n')
  cat('Data Points:',x$data.points,'\n')
  cat('Missing:',x$data.missing,'\n')
  invisible(y)
}

#' cocobot class print method
#' @param x cocobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
#' @method print cocobot
#' @author Charles Dupont
print.cocobot <- print.cobot
