#' cocobot class print method
#' @param x cocobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
#' @method print cocobot
#' @author Charles Dupont

print.cocobot <- function(x, ...) {
  y <- matrix(nrow=length(x$TS),ncol=3)
  dims <- character(length(x$TS))
  for (i in 1:length(x$TS)){
    y[i,] <- c(x$TS[[i]]$ts, sqrt(x$TS[[i]]$var),x$TS[[i]]$pval)
    dims[i] <- x$TS[[i]]$label
  }
  dimnames(y) <- list(dims,c('est','stderr','p'))
  invisible(print(y,...))
}
