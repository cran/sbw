# Plot output from sbwaux
.plot.sbwaux = function(x, plot_cov, ask, ...) {
  if (class(x) != "sbwaux") {
    warning("Object not of class \"sbwaux\"")
    return(invisible(NULL))
  }
  object = x
  weights = object$dat_weights$sbw_weights
  object$dat_weights$sbw_weights = NULL

  par(ask = ask)
  boxplot(weights, main = "Distribution of the weights")

  dat = object$dat_weights
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  bal_tar = object$balance_parameters$bal_tar
  bal_tol = object$balance_parameters$bal_tol
  bal_cov = names(bal_tar)
  if (!is.null(plot_cov)) {
    temp = match(plot_cov, bal_cov)
    check_cov = plot_cov[is.na(temp)]
    if (length(check_cov) > 0) {
      stop(paste(check_cov, "is not found in bal_cov. "))
    } else {
      bal_tar = bal_tar[temp]
      bal_tol = bal_tol[temp]
      bal_cov = bal_cov[temp]
    }
  }
    
  for (i in 1:length(bal_cov)) {
    den_b = density(dat[, bal_cov[i]])
    den_a = spatstat::unnormdensity(dat[, bal_cov[i]], weights = weights)
    max_y = max(den_b$y, den_a$y)
    plot(den_b$x, den_b$y, type = 'l', lwd = 1, lty = 3, col = "gray48",
         ylim = range(c(0, max_y)),
         ylab = 'Density', xlab = bal_cov[i], main = "Balance")
    abline(v = mean(dat[, bal_cov[i]]), lty = 3, col = "gray48")
    lines(den_a$x, den_a$y, col = "gray48", lwd = 1,
          ylab = 'Density', xlab = bal_cov[i])
    abline(v = sum(dat[, bal_cov[i]]*weights), col = "gray48")
    segments(x0 = bal_tar[i], y0 = 0, y1 = max_y/3, col = "black", lty = 4)
    segments(x0 = bal_tar[i] + bal_tol[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    segments(x0 = bal_tar[i] - bal_tol[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    legend("topright", c("Before", "After", "Target", expression(T %+-% tol)),
           col=c("gray48", "gray48", "black", "black"),
           lty = c(3, 1, 4, 2), lwd = 1, cex = 0.75)
    par(ask = ask)
  }
  par(ask = FALSE)
}


# Plot output from sbwcau
.plot.sbwcau = function(x, plot_cov, ask, ...) {
  if (class(x) != "sbwcau") {
    warning("Object not of class \"sbwcau\"")
    return(invisible(NULL))
  }
  object = x
  ind = object$ind
  tre_ind = object$dat_weights[, ind]
  weights0 = object$dat_weights$sbw_weights*(1 - as.numeric(as.character(object$dat_weights[, ind])))
  weights1 = object$dat_weights$sbw_weights*as.numeric(as.character(object$dat_weights[, ind]))

  object$dat_weights$sbw_weights = NULL

  par(ask = ask)
  par(mfrow=c(1, 2))
  boxplot(weights1[tre_ind == 1], main = "Weights in the treated sample")
  boxplot(weights0[tre_ind == 0], main = "Weights in the control sample")
  par(mfrow=c(1, 1))
  
  dat = object$dat_weights
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  
  bal_cov = bal$bal_cov
  
  if (object$par$par_est %in% c("ate", "cate")) {
    bal_tar0 = object$balance_parameters[[1]]$bal_tar
    bal_tol0 = object$balance_parameters[[1]]$bal_tol
    bal_tar1 = object$balance_parameters[[2]]$bal_tar
    bal_tol1 = object$balance_parameters[[2]]$bal_tol
  } else if (object$par$par_est == "att") {
    bal_tar0 = object$balance_parameters$bal_tar
    bal_tol0 = object$balance_parameters$bal_tol
    bal_tar1 = object$balance_parameters$bal_tar
    bal_tol1 = rep(0, length(bal_tar1))
  } else if (object$par$par_est == "atc") {
    bal_tar1 = object$balance_parameters$bal_tar
    bal_tol1 = object$balance_parameters$bal_tol
    bal_tar0 = object$balance_parameters$bal_tar
    bal_tol0 = rep(0, length(bal_tar0))
  }
  
  if (!is.null(plot_cov)) {
    temp = match(plot_cov, bal_cov)
    check_cov = plot_cov[is.na(temp)]
    if (length(check_cov) > 0) {
      stop(paste(check_cov, "is not found in bal_cov. "))
    } else {
      bal_tar0 = bal_tar0[temp]
      bal_tar1 = bal_tar1[temp]
      bal_tol0 = bal_tol0[temp]
      bal_tol1 = bal_tol1[temp]
      bal_cov = bal_cov[temp]
    }
  }
    
  for (i in 1:length(bal_cov)) {
    den_b0 = density(dat[dat[, ind] == 0, bal_cov[i]])
    den_b1 = density(dat[dat[, ind] == 1, bal_cov[i]])
    den_a0 = spatstat::unnormdensity(dat[, bal_cov[i]], weights = weights0)
    den_a1 = spatstat::unnormdensity(dat[, bal_cov[i]], weights = weights1)
    max_y = max(den_b0$y, den_b1$y, den_a0$y, den_a1$y)
    par(mfrow=c(1,2))
    plot(den_b1$x, den_b1$y, type = 'l', lwd = 1, lty = 3,
         ylim = range(c(0, max_y)), col = "gray48", 
         ylab = 'Density', xlab = bal_cov[i], main = "Balance in the treated sample", cex = 0.75)
    abline(v = mean(dat[dat[, ind] == 1, bal_cov[i]]), lwd = 1, lty = 3, col = "gray48")
    lines(den_a1$x, den_a1$y, col = "gray48", lwd = 1,
          ylab = 'Density', xlab = bal_cov[i])
    abline(v = sum(x = dat[, bal_cov[i]]*weights1), lwd = 1, col = "gray48")
    segments(x0 = bal_tar1[i], y0 = 0, y1 = max_y/3, col = "black", lwd = 1, lty = 4)
    segments(x0 = bal_tar1[i] + bal_tol1[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lwd = 1, lty = 2)
    segments(x0 = bal_tar1[i] - bal_tol1[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lwd = 1, lty = 2)

    legend("topright", c("Before", "After", "Target", expression(T%+-% tol)),
           col=c("gray48", "gray48", "black", "black"),
           lty = c(3, 1, 4, 2),
           lwd = c(1, 1, 1, 1), cex = 0.5)
    
    plot(den_b0$x, den_b0$y, type = 'l', lwd = 1, lty = 3,
         ylim = range(c(0, max_y)), 
         ylab = 'Density', xlab = bal_cov[i], main = "Balance in the control sample", cex = 0.75)
    abline(v = mean(dat[dat[, ind] == 0, bal_cov[i]]), lty = 3, col = "gray48")
    lines(den_a0$x, den_a0$y, col = "gray48", lwd = 1,
          ylab = 'Density', xlab = bal_cov[i])
    abline(v = sum(x = dat[, bal_cov[i]]*weights0), col = "gray48")
    segments(x0 = bal_tar0[i], y0 = 0, y1 = max_y/3, col = "black", lty = 4)
    segments(x0 = bal_tar0[i]+bal_tol0[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    segments(x0 = bal_tar0[i]-bal_tol0[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    
    legend("topright", c("Before", "After", "Target", expression(T%+-% tol)),
           col=c("gray48", "gray48", "black", "black"),
           lty = c(3, 1, 4, 2),
           lwd = c(1, 1, 1, 1), cex = 0.5)
    par(ask = ask)
    par(mfrow=c(1, 1))
  }
  par(ask = FALSE)
# }
}

# Plot output from sbwpop
.plot.sbwpop = function(x, plot_cov, ask, ...) {
  if (class(x) != "sbwpop") {
    warning("Object not of class \"sbwpop\"")
    return(invisible(NULL))
  }
  object = x
  ind = object$ind
  tre_ind = object$dat_weights[, ind]
  
  weights0 = object$dat_weights$sbw_weights
  object$dat_weights$sbw_weights = NULL

  par(ask = ask)
  boxplot(weights0[tre_ind == 0], main = "Weights in the complete sample")

  dat = object$dat_weights
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  bal_tar0 = object$balance_parameters$bal_tar
  bal_tol0 = object$balance_parameters$bal_tol

  bal_cov = bal$bal_cov
  
  if (!is.null(plot_cov)) {
    temp = match(plot_cov, bal_cov)
    check_cov = plot_cov[is.na(temp)]
    if (length(check_cov) > 0) {
      stop(paste(check_cov, "is not found in bal_cov. "))
    } else {
      bal_tar0 = bal_tar0[temp]
      bal_tol0 = bal_tol0[temp]
      bal_cov = bal_cov[temp]
      }
    }
    
  for (i in 1:length(bal_cov)) {
    den_b0 = density(dat[dat[, ind] == 0, bal_cov[i]])
    den_a0 = spatstat::unnormdensity(dat[, bal_cov[i]], weights = weights0)
    max_y = max(den_b0$y, den_a0$y)
    plot(den_b0$x, den_b0$y, type = 'l', lwd = 1, lty = 3, col = "gray48",
         ylim = range(c(0, max_y)),
         ylab = 'Density', xlab = bal_cov[i], main = "Balance in the complete sample")
    abline(v = mean(dat[dat[, ind] == 0, bal_cov[i]]), lty = 3, col = "gray48")
    lines(den_a0$x, den_a0$y, col = "gray48", lwd = 1,
          ylab = 'Density', xlab = bal_cov[i])
    abline(v = sum(x = dat[, bal_cov[i]]*weights0), col = "gray48")
    segments(x0 = bal_tar0[i], y0 = 0, y1 = max_y/3, col = "black", lty = 4)
    segments(x0 = bal_tar0[i] + bal_tol0[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    segments(x0 = bal_tar0[i] - bal_tol0[i], y0 = 2*max_y/3, y1 = max_y, col = "black", lty = 2)
    legend("topright", c("Before", "After", "Target", expression(T %+-% tol)),
           col=c("gray48", "gray48", "black", "black"),
           lty = c(3, 1, 4, 2),
           lwd = c(1, 1, 1, 1), cex = 0.75)
    par(ask = ask)
  }
  par(ask = FALSE)
}

#' Visualize output from \code{sbw}
#'
#' @description Function for visualizing the output from \code{\link[sbw]{sbw}}.
#'
#' @param object an object from function \code{\link[sbw]{sbw}}.
#' @param plot_cov names of covariates for which balance is to be displayed.  If \code{NULL}, all of the covariates will be displayed.
#' @param ask logical. If \code{TRUE} (and the R session is interactive) the user is asked for input, before a new figure is drawn.
#' @param ... ignored arguments.
#' 
#' @importFrom spatstat unnormdensity
#' 
#' @examples 
#' # Please see the examples in the function sbw above.
#' @export
#' 
visualize = function(object, plot_cov, ask = TRUE, ...) {
  if (missing(plot_cov)) plot_cov = NULL
  if (class(object) == "sbwaux") {
    .plot.sbwaux(x = object, plot_cov = plot_cov, ask = ask, ...)
  } else if (class(object) == "sbwcau") {
    .plot.sbwcau(x = object, plot_cov = plot_cov, ask = ask, ...)
  } else if (class(object) == "sbwpop") {
    .plot.sbwpop(x = object, plot_cov = plot_cov, ask = ask, ...)
  } else stop("Please use one of the calls from sbw.")
}
