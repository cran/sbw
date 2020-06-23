# Summarize output from sbwaux
.summary.sbwaux = function(object, digits, ...) {
  if (class(object) != "sbwaux") {
    warning("Object not of class \"sbwaux\".")
    return(invisible(NULL))
  }

  weights = object$dat_weights$sbw_weights
  object$dat_weights$sbw_weights = NULL
  var_weights = var(weights)
  cv_weights = sd(weights)/mean(weights)
  effective_sample_size = object$effective_sample_size
  dat = object$dat_weights[,object$bal$bal_cov]
  bal_tar = object$balance_parameters$bal_tar
  bal_tol = object$balance_parameters$bal_tol
  # bal_cov = names(bal_tar)
  # bal_std = object$balance_parameters$bal_std

  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  dat = data.matrix(dat)
  storage.mode(dat) = "numeric"
  colnames(dat) = object$bal$bal_cov
  
  means_b = colMeans(dat)
  sds_b = apply(dat, 2, sd)

  means_a = as.vector(t(weights)%*%as.matrix(dat))
  target = rep(NA, length(means_b))
  target[match(names(bal_tar),colnames(dat))] = bal_tar
  tolerance = rep(NA, length(means_b))
  if (object$bal$bal_std %in% c("group", "manual")) {
    tolerance[match(names(bal_tar),colnames(dat))] = object$balance_parameters$bal_tol_ori
  } else if (object$bal$bal_std %in% "target") {
    tolerance[match(names(bal_tar),colnames(dat))] = object$balance_parameters$bal_tol
  }
  
  if (object$bal$bal_std %in% "group") {
    dif_b = as.vector(abs(target - means_b)/sds_b)
  } else if (object$bal$bal_std %in% "manual") {
    dif_b = as.vector(abs(target - means_b))
  } else if (object$bal$bal_std %in% "target") {
    dif_b = as.vector(abs(target - means_b))
  }

  if (object$bal$bal_std %in% "group") {
    dif_a = as.vector(abs(target - means_a)/sds_b)
  } else if (object$bal$bal_std %in% "manual") {
    dif_a = as.vector(abs(target - means_a))
  } else if (object$bal$bal_std %in% "target") {
    dif_a = as.vector(abs(target - means_a))
  }
  
  tab = cbind(paste(round(means_b, digits)," / ", round(dif_b, digits), sep=""),
              paste(round(means_a, digits)," / ", round(dif_a, digits), sep=""),
              paste(round(target, digits)), paste(round(tolerance, digits)))

  rownames(tab) = names(means_b)
  colnames(tab) =  c("Before", "After", "Target", "Tolerance")
  
  tab_1 = round(cbind(means_b, means_a, target), digits = digits)
  rownames(tab_1) = names(means_b)
  colnames(tab_1) = c("Before", "After", "Target")
  
  tab_2 = round(cbind(dif_b, dif_a, tolerance), digits = digits)
  rownames(tab_2) = names(means_b)
  colnames(tab_2) = c("Before", "After", "Tolerance")
  
  if (!is.null(object$shadow_price)) {
    colnames(object$shadow_price) = c("Upper", "Lower")
  }
  
  cat("\n")
  cat("Variance of the weights in the sample: ", format(var_weights, digits = digits), "\n")
  cat("\n")
  cat("Coefficient of variation of the weights in the sample: ", format(cv_weights, digits = digits),"\n")
  cat("\n")
  cat("Effective sample size: ", format(effective_sample_size, digits = digits),"\n")
  cat("\n")
  cat("Means of the sample before and after weighting: ", "\n")
  print(tab_1, digits = digits)
  cat("\n")
  cat("TASDM of the sample before and after weighting: ", "\n")
  print(tab_2, digits = digits)
  cat("\n")
  cat("Shadow prices for the sample: ", "\n")
  print(object$shadow_price, digits = digits)
  cat("\n")
  
  invisible(list(variance = var_weights, coefficient_variation = cv_weights, 
                 effective_sample_size = effective_sample_size,
                 balance_table = list(mean = tab_1, TASDM = tab_2),
                 shadow_price = object$shadow_price))
}


# Summarize output from sbwcau
.summary.sbwcau = function(object, digits, ...) {
  if (class(object) != "sbwcau") {
    warning("Object not of class \"sbwcau\"")
    return(invisible(NULL))
  }
  ind = object$ind
  out = object$out
  tre_ind = as.numeric(as.character(object$dat_weights[, ind]))
  weights0 = object$dat_weights$sbw_weights*(1 - tre_ind)
  weights1 = object$dat_weights$sbw_weights*tre_ind

  object$dat_weights$sbw_weights = NULL
  var_weights0 = var(weights0[tre_ind == 0])
  var_weights1 = var(weights1[tre_ind == 1])
  cv_weights0 = sd(weights0[tre_ind == 0])/mean(weights0[tre_ind == 0])
  cv_weights1 = sd(weights1[tre_ind == 1])/mean(weights1[tre_ind == 1])
  effective_sample_size0 = sum(weights0[tre_ind == 0])^2/sum(weights0[tre_ind == 0]^2)
  effective_sample_size1 = sum(weights1[tre_ind == 1])^2/sum(weights1[tre_ind == 1]^2)
  
  variance = c(var_weights1, var_weights0)
  names(variance) = c("treated", "control")
  coef_var = c(cv_weights1, cv_weights0)
  names(coef_var) = c("treated", "control")
  effective_sample_size = c(effective_sample_size1, effective_sample_size0)
  names(effective_sample_size) = c("treated", "control")
  
  dat = object$dat_weights
  dat = object$dat_weights[,object$bal$bal_cov]
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  dat = data.matrix(dat)
  storage.mode(dat) = "numeric"
  colnames(dat) = object$bal$bal_cov
    
  # bal_cov = bal$bal_cov
  means_b0 = colMeans(dat[tre_ind == 0,,drop=FALSE])
  sds_b0 = apply(dat[tre_ind == 0,,drop=FALSE], 2, sd)
  means_b1 = colMeans(dat[tre_ind == 1,,drop=FALSE])
  sds_b1 = apply(dat[tre_ind == 1,,drop=FALSE], 2, sd)
  sds_b = apply(dat, 2, sd)

  means_a0 = as.vector(t(weights0)%*%as.matrix(dat))
  means_a1 = as.vector(t(weights1)%*%as.matrix(dat))

  target = rep(NA, length(means_b0))
  tolerance = rep(NA, length(means_b0))
  
  target[match(bal$bal_cov,colnames(dat))] = object$bal$bal_tar
  tolerance[match(bal$bal_cov,colnames(dat))] = object$bal$bal_tol
  
  if (object$bal$bal_std %in% "group") {
    dif_b0 = as.vector(abs(target - means_b0)/sds_b0)
    dif_b1 = as.vector(abs(target - means_b1)/sds_b1)
  } else if (object$bal$bal_std %in% "target") {
    if (object$par$par_est %in% c("ate", "cate")) {
      dif_b0 = as.vector(abs(target - means_b0)/sds_b)
      dif_b1 = as.vector(abs(target - means_b1)/sds_b)
    } else if (object$par$par_est %in% c("att", "atc")) {
      dif_b0 = as.vector(abs(target - means_b0)/sds_b1)
      dif_b1 = as.vector(abs(target - means_b1)/sds_b0)
    }
  } else if (object$bal$bal_std %in% "manual") {
    dif_b0 = as.vector(abs(target - means_b0))
    dif_b1 = as.vector(abs(target - means_b1))
  }

  tab_b = cbind(paste(round(means_b1, digits), " / ", round(dif_b1, digits), sep=""),
                paste(round(means_b0, digits), " / ", round(dif_b0, digits), sep=""),
                paste(round(target, digits)), paste(round(tolerance, digits)))
  rownames(tab_b) = names(means_b0)
  colnames(tab_b) =  c("Treated", "Control", "Target", "Tolerance")

  tab_b1 = round(cbind(means_b1, means_b0, target), digits = digits)
  rownames(tab_b1) = names(means_b0)
  colnames(tab_b1) = c("Treated", "Control", "Target")

  tab_b2 = round(cbind(dif_b1, dif_b0, tolerance), digits = digits)
  rownames(tab_b2) = names(means_b0)
  colnames(tab_b2) = c("Treated", "Control", "Tolerance")
  
  if (object$bal$bal_std %in% "group") {
    dif_a0 = as.vector(abs(target - means_a0)/sds_b0)
    dif_a1 = as.vector(abs(target - means_a1)/sds_b1)
  } else if (object$bal$bal_std %in% "target") {
    if (object$par$par_est %in% c("ate", "cate")) {
      dif_a0 = as.vector(abs(target - means_a0)/sds_b)
      dif_a1 = as.vector(abs(target - means_a1)/sds_b)
    } else if (object$par$par_est %in% c("att", "atc")) {
      dif_a0 = as.vector(abs(target - means_a0)/sds_b1)
      dif_a1 = as.vector(abs(target - means_a1)/sds_b0)
    }
  } else if (object$bal$bal_std %in% "manual") {
    dif_a0 = as.vector(abs(target - means_a0))
    dif_a1 = as.vector(abs(target - means_a1))
  }
  
  tab_a = cbind(paste(round(means_a1, digits)," / ", round(dif_a1, digits), sep=""),
                paste(round(means_a0, digits)," / ", round(dif_a0, digits), sep=""),
                paste(round(target, digits)), paste(round(tolerance, digits)))
  rownames(tab_a) = names(means_b0)
  colnames(tab_a) = c("Treated", "Control", "Target", "Tolerance")

  tab_a1 = round(cbind(means_a1, means_a0, target), digits = digits)
  rownames(tab_a1) = names(means_b0)
  colnames(tab_a1) = c("Treated", "Control", "Target")

  tab_a2 = round(cbind(dif_a1, dif_a0, tolerance), digits = digits)
  rownames(tab_a2) = names(means_b0)
  colnames(tab_a2) = c("Treated", "Control", "Tolerance")
  
  cat("\n")
  cat("Variance of the weights in the treated sample: ", format(var_weights1, digits = digits), "\n")
  cat("Variance of the weights in the control sample: ", format(var_weights0, digits = digits), "\n")
  cat("\n")
  cat("Coefficient of variation of the weights in the treated sample: ", format(cv_weights1, digits = digits),"\n")
  cat("Coefficient of variation of the weights in the control sample: ", format(cv_weights0, digits = digits),"\n")
  cat("\n")
  cat("Effective sample size of the treated sample: ", format(effective_sample_size1, digits = digits),"\n")
  cat("Effective sample size of the control sample: ", format(effective_sample_size0, digits = digits),"\n")
  cat("\n")
  cat("Means of the treated and control samples before weighting: ", "\n")
  print(tab_b1, digits = digits)
  cat("\n")
  cat("TASDM of the treated and control samples before weighting: ", "\n")
  print(tab_b2, digits = digits)
  cat("\n")
  cat("Means of the treated and control samples after weighting: ", "\n")
  print(tab_a1, digits = digits)
  cat("\n")
  cat("TASDM of the treated and control samples after weighting: ", "\n")
  print(tab_a2, digits = digits)
  cat("\n")
  if (object$par$par_est %in% c("ate", "cate")) {
    cat("Shadow prices for the treated sample: ", "\n")
    print(object$shadow_price[[2]], digits = digits)
    cat("\n")
    cat("Shadow prices for the control sample: ", "\n")
    print(object$shadow_price[[1]], digits = digits)
    cat("\n")
  } else if (object$par$par_est == "att") {
    cat("Shadow prices for the control sample: ", "\n")
    print(object$shadow_price, digits = digits)
    cat("\n")
  } else if (object$par$par_est == "atc") {
    cat("Shadow prices for the treated sample: ", "\n")
    print(object$shadow_price, digits = digits)
    cat("\n")
  }
  
  invisible(list(variance = variance, coefficient_variation = coef_var, effective_sample_size = effective_sample_size,
                 balance_table = list(mean_before = tab_b1, TASDM_before = tab_b2,
                                      mean_after = tab_a1, TASDM_after = tab_a2), 
                 shadow_price = object$shadow_price))
}


# Summarize output from sbwpop
.summary.sbwpop = function(object, digits, ...) {
  if (class(object) != "sbwpop") {
    warning("Object not of class \"sbwpop\"")
    return(invisible(NULL))
  }
  ind = object$ind
  out = object$out
  tre_ind = as.numeric(as.character(object$dat_weights[, ind]))
  
  weights0 = object$dat_weights$sbw_weights
  object$dat_weights$sbw_weights = NULL
  var_weights0 = var(weights0[tre_ind == 0])
  cv_weights0 = sd(weights0[tre_ind == 0])/mean(weights0[tre_ind == 0])
  effective_sample_size0 = sum(weights0[tre_ind == 0])^2/sum(weights0[tre_ind == 0]^2)
  
  dat = object$dat_weights
  dat = object$dat_weights[,object$bal$bal_cov]
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  dat = data.matrix(dat)
  storage.mode(dat) = "numeric"
  colnames(dat) = object$bal$bal_cov
  
  #bal_cov = bal$bal_cov
  means_b0 = colMeans(dat[tre_ind == 0,,drop=FALSE])
  sds_b0 = apply(dat[tre_ind == 0,,drop=FALSE], 2, sd)
  means_b1 = colMeans(dat)
  sds_b1 = apply(dat, 2, sd)

  temp = dat
  temp[is.na(temp)] = 0
  means_a0 = as.vector(t(weights0)%*%as.matrix(temp))
  means_a1 = means_b1
  sds_a1 = sds_b1

  target = rep(NA, length(means_b0))
  tolerance = rep(NA, length(means_b0))
  
  target[match(bal$bal_cov,colnames(dat))] = object$bal$bal_tar
  tolerance[match(bal$bal_cov,colnames(dat))] = object$bal$bal_tol

  if (object$bal$bal_std %in% "group") {
    dif_b0 = as.vector(abs(target - means_b0)/sds_b0)
    dif_b1 = as.vector(abs(target - means_b1)/sds_b1)
  } else if (object$bal$bal_std %in% "target") {
    dif_b0 = as.vector(abs(target - means_b0)/sds_b1)
    dif_b1 = as.vector(abs(target - means_b1)/sds_b1)
  } else if (object$bal$bal_std %in% "manual") {
    dif_b0 = as.vector(abs(target - means_b0))
    dif_b1 = as.vector(abs(target - means_b1))
  }
  
  tab_b = cbind(paste(round(means_b0, digits)," / ", round(dif_b0, digits), sep=""),
                paste(round(means_b1, digits)," / ", round(dif_b1, digits), sep=""),
                paste(round(target, digits)), paste(round(tolerance, digits)))
  rownames(tab_b) = names(means_b0)
  colnames(tab_b) =  c("Non-missing", "Total", "Target", "Tolerance")
  
  tab_b1 = round(cbind(means_b0, means_b1, target), digits = digits)
  rownames(tab_b1) = names(means_b0)
  colnames(tab_b1) = c("Non-missing", "Total", "Target")
  
  tab_b2 = round(cbind(dif_b0, dif_b1, tolerance), digits = digits)
  rownames(tab_b2) = names(means_b0)
  colnames(tab_b2) = c("Non-missing", "Total", "Tolerance")
  
  if (object$bal$bal_std %in% "group") {
    dif_a0 = as.vector(abs(target - means_a0)/sds_b0)
    dif_a1 = as.vector(abs(target - means_a1)/sds_b1)
  } else if (object$bal$bal_std %in% "target") {
    dif_a0 = as.vector(abs(target - means_a0)/sds_b1)
    dif_a1 = as.vector(abs(target - means_a1)/sds_b1)
  } else if (object$bal$bal_std %in% "manual") {
    dif_a0 = as.vector(abs(target - means_a0))
    dif_a1 = as.vector(abs(target - means_a1))
  }
  
  tab_a = cbind(paste(round(means_a0, digits), " / ", round(dif_a0, digits), sep=""),
                paste(round(means_a1, digits), " / ", round(dif_a1, digits), sep=""),
                paste(round(target, digits)), paste(round(tolerance, digits)))
  rownames(tab_a) = names(means_b0)
  colnames(tab_a) =  c("Non-missing", "Total", "Target", "Tolerance")

  tab_a1 = round(cbind(means_a0, means_a1, target), digits = digits)
  rownames(tab_a1) = names(means_b0)
  colnames(tab_a1) = c("Non-missing", "Total", "Target")
  
  tab_a2 = round(cbind(dif_a0, dif_a1, tolerance), digits = digits)
  rownames(tab_a2) = names(means_b0)
  colnames(tab_a2) = c("Non-missing", "Total", "Tolerance")

  cat("\n")
  cat("Variance of the weights in the non-missing sample: ", format(var_weights0, digits = digits), "\n")
  cat("\n")
  cat("Coefficient of variation of the weights in the non-missing sample: ", format(cv_weights0, digits = digits),"\n")
  cat("\n")
  cat("Effective sample size of the non-missing sample: ", format(effective_sample_size0, digits = digits),"\n")
  cat("\n")
  cat("Means of the non-missing and total sample before weighting: ", "\n")
  print(tab_b1, digits = digits)
  cat("\n")
  cat("TASDM of the non-missing and total sample before weighting: ", "\n")
  print(tab_b2, digits = digits)
  cat("\n")
  cat("Means of the non-missing and total sample after weighting: ", "\n")
  print(tab_a1, digits = digits)
  cat("\n")
  cat("TASDM of the non-missing and total sample after weighting: ", "\n")
  print(tab_a2, digits = digits)
  cat("\n")
  cat("Shadow prices for the non-missing sample: ", "\n")
  print(object$shadow_price, digits = digits)
  cat("\n")

  invisible(list(variance = var_weights0, coefficient_variation = cv_weights0, 
                 balance_table = list(mean_before = tab_b1, TASDM_before = tab_b2,
                                      mean_after = tab_a1, TASDM_after = tab_a2), 
                 shadow_price = object$shadow_price))
}

#' Summarize output from \code{sbw}
#'
#' @description Function for summarizing the output from \code{\link[sbw]{sbw}}.
#'
#' @param object an object from the class \code{sbwcau} or \code{sbwpop} obtained after using \code{\link[sbw]{sbw}}.
#' @param digits The number of significant digits that will be displayed. The default is \code{6}.
#' @param ... ignored arguments.
#' 
#' @importFrom spatstat unnormdensity
#' 
#' @return A list with the following elements:
#' @return \code{variance}{, variance of the weights}
#' @return \code{coefficient_variation}{, coefficient of variation of the weights}
#' @return \code{effective_sample_size}{, effective sample size}
#' @return \code{balance_table}{, mean/TASDM balance tables for samples before/after weighting}
#' @return \code{shadow_price}{, dual tables or shadow prices for the balanced groups}
#' 
#' @examples 
#' # Please see the examples in the function sbw above.
#' @export
#' 
summarize = function(object, digits = 6, ...) {
  if (class(object) == "sbwaux") {
    .summary.sbwaux(object, digits = digits, ...)
  } else if (class(object) == "sbwcau") {
    .summary.sbwcau(object, digits = digits, ...)
  } else if (class(object) == "sbwpop") {
    .summary.sbwpop(object, digits = digits, ...)
  } else stop("Please use one of the calls from sbw.")
}