# Estimate output from sbwcau
.estimate.sbwcau = function(object, digits = 2, ...) {
  if (class(object) != "sbwcau") {
    warning("Object not of class \"sbwcau\"")
    return(invisible(NULL))
  }
  ind = object$ind
  out = object$out
  # if (is.null(out)) {stop("argument \"out\" is missing in the function \"sbw\".")}
  dat = object$dat_weights
  if (sum(1 - is.na(match(out, colnames(dat)))) == 0) {
    stop("Please specify a correct string for out.")
  }
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  
  if (object$par$par_est == "att") {
    tre_ind = dat[, ind]
    weights0 = dat$weights*(1 - tre_ind)
    weights1 = dat$weights*tre_ind
    n = sum(tre_ind == 1)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
     
    var_cau = colMeans(as.matrix((as.matrix(n*weights1*Y - sum(weights1*Y) 
                        - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))[tre_ind == 1,])^2) 
    + colMeans(as.matrix((as.matrix(n*weights0*Y - sum(weights0*Y)
                        - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))[tre_ind == 1,])^2)
    sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "atc") {
    tre_ind = dat[, ind]
    weights0 = dat$weights*(1 - tre_ind)
    weights1 = dat$weights*tre_ind
    n = sum(tre_ind == 0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans(as.matrix((as.matrix(n*weights1*Y - sum(weights1*Y) 
                        - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))[tre_ind == 0,])^2)
    +  colMeans(as.matrix((as.matrix(n*weights0*Y - sum(weights0*Y) 
                        - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))[tre_ind == 0,])^2)
    sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "ate") {
    tre_ind = dat[, ind]
    weights0 = dat$weights*(1 - tre_ind)
    weights1 = dat$weights*tre_ind
    n = length(weights0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$weights = NULL
    dat = dat[, object$bal$bal_cov]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans((as.matrix(n*weights1*Y - sum(weights1*Y) 
                   - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))^2)
            + colMeans((as.matrix(n*weights0*Y - sum(weights0*Y) 
                   - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))^2)
    sd_cau = sqrt(var_cau/n)
  }
  if (object$par$par_est == "cate") {
    if (class(object$par$par_tar) == "character") {
      dat = subset(dat, eval(parse(text = object$par$par_tar)))
      dat[fac_ind] = NULL
    } else if (class(object$par$par_tar) == "numeric") {
      if (sum(fac_ind) >= 1) {
        dat = dat[apply(dat[fac_ind] == object$par$par_tar[match(names(fac_ind), names(object$par$par_tar))][fac_ind], 1, prod, na.rm = TRUE) %in% 1,]
      }
      dat[fac_ind] = NULL
    }
    tre_ind = dat[, ind]
    weights0 = dat$weights*(1 - tre_ind)
    weights1 = dat$weights*tre_ind
    n = length(weights0)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$weights = NULL
    dat = dat[, object$bal$bal_cov[!(object$bal$bal_cov %in% names(which(fac_ind == TRUE)))]]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_cau = colMeans((as.matrix(n*weights1*Y - sum(weights1*Y) 
                    - dat%*%solve(t(weights1*dat)%*%dat)%*%(t(weights1*dat)%*%Y)*(n*weights1 - 1)))^2)
            + colMeans((as.matrix(n*weights0*Y - sum(weights0*Y) 
                    - dat%*%solve(t(weights0*dat)%*%dat)%*%(t(weights0*dat)%*%Y)*(n*weights0 - 1)))^2)
    sd_cau = sqrt(var_cau/n)
  }
  
  estimates = crossprod(weights1 - weights0, Y)
  estimates = as.vector(estimates)
  
  cau_table = cbind(estimates, sd_cau, estimates/sd_cau, 
                    pt(q = abs(estimates/sd_cau), df = n - 1, lower.tail = FALSE),
                    estimates + sd_cau*qt(0.025, n-1),
                    estimates + sd_cau*qt(0.975, n-1))
  rownames(cau_table) = out
  colnames(cau_table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "2.5 %", "97.5 %")
  cat("\n")
  cat(paste(toupper(object$par$par_est), ": ", sep = ""), "\n")
  print(cau_table, digits = digits)
  cat("\n")
  invisible(list(cau_table = cau_table))
}

# Estimate output from sbwpop
.estimate.sbwpop = function(object, digits = 2, ...) {
  if (class(object) != "sbwpop") {
    warning("Object not of class \"sbwpop\"")
    return(invisible(NULL))
  }
  if (object$par$par_est == "pop") {
    ind = object$ind
    out = object$out
    # if (is.null(out)) {stop("argument \"out\" is missing in the function \"sbw\".")}
    dat = object$dat_weights
    if (sum(1 - is.na(match(out, colnames(dat)))) == 0) {
      stop("Please specify a correct string for out.")
    }
    fac_ind = sapply(dat, is.factor)
    dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
    if (class(object$par$par_tar) == "character") {
      dat = subset(dat, eval(parse(text = object$par$par_tar)))
      dat[fac_ind] = NULL
    } else if (class(object$par$par_tar) == "numeric") {
      if (sum(fac_ind) >= 1) {
        dat = dat[apply(dat[fac_ind] == object$par$par_tar[match(names(fac_ind), names(object$par$par_tar))][fac_ind], 1, prod, na.rm = TRUE) %in% 1,]
      }
      dat[fac_ind] = NULL
    }
    tre_ind = dat[, ind]
    dat = dat[tre_ind == 0,]
    weights = dat$weights
    n = length(weights)
    Y = as.matrix(dat[, out])
    dat[, ind] = NULL
    dat[, out] = NULL
    dat$weights = NULL
    dat = dat[, object$bal$bal_cov[!(object$bal$bal_cov %in% names(which(fac_ind == TRUE)))]]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_pop = colMeans((as.matrix(n*weights*Y - sum(weights*Y) 
                    - dat%*%solve(t(weights*dat)%*%dat)%*%(t(weights*dat)%*%Y)*(n*weights - 1)))^2)
    sd_pop = sqrt(var_pop/n)
  }
  
  estimates = crossprod(weights, Y)
  estimates = as.vector(estimates)
  pop_table = cbind(estimates, sd_pop, estimates/sd_pop, 
                    pt(q = abs(estimates/sd_pop), df = n - 1, lower.tail = FALSE),
                    estimates + sd_pop*qt(0.025, n-1),
                    estimates + sd_pop*qt(0.975, n-1))
  rownames(pop_table) = out
  colnames(pop_table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "2.5 %", "97.5 %")
  cat("\n")
  cat(paste(toupper(object$par$par_est), ": ", sep = ""), "\n")
  print(pop_table, digits = digits)
  cat("\n")
  invisible(list(pop_table = pop_table))
}

# Estimate output from sbwaux
.estimate.sbwaux = function(object, digits = 2, ...) {
  if (class(object) != "sbwaux") {
    warning("Object not of class \"sbwaux\"")
    return(invisible(NULL))
  }
  if (object$par$par_est == "aux") {
    out = object$out
    # if (is.null(out)) {stop("argument \"out\" is missing in the function \"sbw\".")}
    dat = object$dat_weights
    if (sum(1 - is.na(match(out, colnames(dat)))) == 0) {
      stop("Please specify a correct string for out.")
    }
    fac_ind = sapply(dat, is.factor)
    dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
    if (class(object$par$par_tar) == "numeric") {
      if (sum(fac_ind) >= 1) {
        dat = dat[apply(dat[fac_ind] == object$par$par_tar[match(names(fac_ind), names(object$par$par_tar))][fac_ind], 1, prod, na.rm = TRUE) %in% 1,]
      }
      dat[fac_ind] = NULL
    }
    weights = dat$weights
    n = length(weights)
    Y = as.matrix(dat[, out])
    dat = dat[, object$bal$bal_cov[!(object$bal$bal_cov %in% names(which(fac_ind == TRUE)))]]
    dat = as.matrix(dat)
    tmp = cor(dat)
    tmp[upper.tri(tmp)] = 0
    diag(tmp) = 0
    dat = dat[,!apply(tmp,2,function(x) any(x > 0.9999))]
    
    var_aux = colMeans((as.matrix(n*weights*Y - sum(weights*Y) 
                                  - dat%*%solve(t(weights*dat)%*%dat)%*%(t(weights*dat)%*%Y)*(n*weights - 1)))^2)
    sd_aux = sqrt(var_aux/n)
  }
  
  estimates = crossprod(weights, Y)
  estimates = as.vector(estimates)
  aux_table = cbind(estimates, sd_aux, estimates/sd_aux, 
                    pt(q = abs(estimates/sd_aux), df = n - 1, lower.tail = FALSE),
                    estimates + sd_aux*qt(0.025, n-1),
                    estimates + sd_aux*qt(0.975, n-1))
  rownames(aux_table) = out
  colnames(aux_table) = c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "2.5 %", "97.5 %")
  cat("\n")
  cat(paste(toupper(object$par$par_est), ": ", sep = ""), "\n")
  print(aux_table, digits = digits)
  cat("\n")
  invisible(list(aux_table = aux_table))
}

#' Estimate causal contrasts and population means
#'
#' @description Function for estimating causal contrasts and population means using the output from \code{\link[sbw]{sbw}}.
#'
#' @param object an object from function \code{\link[sbw]{sbw}}.
#' @param digits a scalar with the number of significant digits used to display the estimates. The default is 6.
#' @param ... ignored arguments
#' 
#' @examples 
#' # Please see the examples in sbw.
#' @export
#' 
estimate = function(object, digits = 6, ...) {
  if (class(object) == "sbwcau") {
    .estimate.sbwcau(object, digits = digits, ...)
  } else if (class(object) == "sbwpop") {
    .estimate.sbwpop(object, digits = digits, ...)
  } else if (class(object) == "sbwaux") {
    .estimate.sbwaux(object, digits = digits, ...)
  } else stop("Please use one of the calls from sbw.")
}
