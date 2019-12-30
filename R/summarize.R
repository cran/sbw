# Summarize output from sbwaux
.summary.sbwaux = function(object, digits = 2, ...) {
  if (class(object) != "sbwaux") {
    warning("Object not of class \"sbwaux\".")
    return(invisible(NULL))
  }

  weights = object$dat_weights$weights
  object$dat_weights$weights = NULL
  var_weights = var(weights)
  cv_weights = sd(weights)/mean(weights)
  eff_size = object$eff_size
  dat = object$dat_weights
  bal_tar = object$target$bal_tar
  bal_tol = object$target$bal_tol
  bal_cov = names(bal_tar)
  # bal_std = object$target$bal_std

  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  means_b = colMeans(dat)
  sds_b = apply(dat, 2, sd)

  means_a = as.vector(t(weights)%*%as.matrix(dat))
  target = rep(NA, length(means_b))
  target[match(names(bal_tar),colnames(dat))] = bal_tar
  tolerance = rep(NA, length(means_b))
  tolerance[match(names(bal_tar),colnames(dat))] = object$target$bal_tol_ori
  
  if (object$bal$bal_std == TRUE) {
    dif_b = as.vector(abs(target - means_b)/sds_b)
  } else if (object$bal$bal_std == FALSE) {
    dif_b = as.vector(abs(target - means_b))
  }

  if (object$bal$bal_std == TRUE) {
    dif_a = as.vector(abs(target - means_a)/sds_b)
  } else if (object$bal$bal_std == FALSE) {
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
  
  colnames(object$dual_table) = c("Upper", "Lower")
  
  cat("\n")
  cat("Variance of the weights in the sample: ", format(var_weights, digits = digits), "\n")
  cat("\n")
  cat("Coefficient of variation of the weights in the sample: ", format(cv_weights, digits = digits),"\n")
  cat("\n")
  cat("Effective sample size: ", format(eff_size, digits = digits),"\n")
  cat("\n")
  cat("Means of the sample before and after weighting: ", "\n")
  print(tab_1, digits = digits)
  cat("\n")
  cat("ASDM of the sample before and after weighting: ", "\n")
  print(tab_2, digits = digits)
  cat("\n")
  cat("Shadow prices for the sample: ", "\n")
  print(object$dual_table, digits = digits)
  cat("\n")
  
  invisible(list(variance = var_weights, coef_variance = cv_weights, 
                 eff_size = eff_size,
                 mean_table = list(means = tab_1, ASDM = tab_2),
                 dual_table = object$dual_table))
}


# Summarize output from sbwcau
.summary.sbwcau = function(object, digits = 2, ...) {
  if (class(object) != "sbwcau") {
    warning("Object not of class \"sbwcau\"")
    return(invisible(NULL))
  }
  ind = object$ind
  out = object$out
  tre_ind = as.numeric(as.character(object$dat_weights[, ind]))
  weights0 = object$dat_weights$weights*(1 - tre_ind)
  weights1 = object$dat_weights$weights*tre_ind

  object$dat_weights$weights = NULL
  var_weights0 = var(weights0[tre_ind == 0])
  var_weights1 = var(weights1[tre_ind == 1])
  cv_weights0 = sd(weights0[tre_ind == 0])/mean(weights0[tre_ind == 0])
  cv_weights1 = sd(weights1[tre_ind == 1])/mean(weights1[tre_ind == 1])
  eff_size0 = sum(weights0[tre_ind == 0])^2/sum(weights0[tre_ind == 0]^2)
  eff_size1 = sum(weights1[tre_ind == 1])^2/sum(weights1[tre_ind == 1]^2)
  
  variance = c(var_weights1, var_weights0)
  names(variance) = c("treated", "control")
  coef_var = c(cv_weights1, cv_weights0)
  names(coef_var) = c("treated", "control")
  eff_size = c(eff_size1, eff_size0)
  names(eff_size) = c("treated", "control")
  
  dat = object$dat_weights
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  dat[, ind] = NULL
  dat[, out] = NULL
  
  bal_cov = bal$bal_cov
  means_b0 = colMeans(dat[tre_ind == 0,])
  sds_b0 = apply(dat[tre_ind == 0,], 2, sd)
  means_b1 = colMeans(dat[tre_ind == 1,])
  sds_b1 = apply(dat[tre_ind == 1,], 2, sd)

  means_a0 = as.vector(t(weights0)%*%as.matrix(dat))
  means_a1 = as.vector(t(weights1)%*%as.matrix(dat))

  target = rep(NA, length(means_b0))
  tolerance = rep(NA, length(means_b0))
  
  target[match(bal_cov,colnames(dat))] = object$bal$bal_tar
  tolerance[match(bal_cov,colnames(dat))] = object$bal$bal_tol
  
  if (object$bal$bal_std == TRUE) {
    dif_b0 = as.vector(abs(target - means_b0)/sds_b0)
    dif_b1 = as.vector(abs(target - means_b1)/sds_b1)
  } else if (object$bal$bal_std == FALSE) {
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
  
  if (object$bal$bal_std == TRUE) {
    dif_a0 = as.vector(abs(target - means_a0)/sds_b0)
    dif_a1 = as.vector(abs(target - means_a1)/sds_b1)
  } else if (object$bal$bal_std == FALSE) {
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
  cat("Effective sample size of the treated sample: ", format(eff_size1, digits = digits),"\n")
  cat("Effective sample size of the control sample: ", format(eff_size0, digits = digits),"\n")
  cat("\n")
  cat("Means of the treated and control samples before weighting: ", "\n")
  print(tab_b1, digits = digits)
  cat("\n")
  cat("ASDM of the treated and control samples before weighting: ", "\n")
  print(tab_b2, digits = digits)
  cat("\n")
  cat("Means of the treated and control samples after weighting: ", "\n")
  print(tab_a1, digits = digits)
  cat("\n")
  cat("ASDM of the treated and control samples after weighting: ", "\n")
  print(tab_a2, digits = digits)
  cat("\n")
  if (object$par$par_est %in% c("ate", "cate")) {
    cat("Shadow prices for the treated sample: ", "\n")
    print(object$dual_table[[2]], digits = digits)
    cat("\n")
    cat("Shadow prices for the control sample: ", "\n")
    print(object$dual_table[[1]], digits = digits)
    cat("\n")
  } else if (object$par$par_est == "att") {
    cat("Shadow prices for the control sample: ", "\n")
    print(object$dual_table, digits = digits)
    cat("\n")
  } else if (object$par$par_est == "atc") {
    cat("Shadow prices for the treated sample: ", "\n")
    print(object$dual_table, digits = digits)
    cat("\n")
  }
  
  invisible(list(variance = variance, coef_variance = coef_var, eff_size = eff_size,
                 mean_table = list(means_before = tab_b1, ASDM_before = tab_b2,
                                   means_after = tab_a1, ASDM_after = tab_a2), 
                 dual_table = object$dual_table))
}


# Summarize output from sbwpop
.summary.sbwpop = function(object, digits = 2, ...) {
  if (class(object) != "sbwpop") {
    warning("Object not of class \"sbwpop\"")
    return(invisible(NULL))
  }
  ind = object$ind
  out = object$out
  tre_ind = as.numeric(as.character(object$dat_weights[, ind]))
  
  weights0 = object$dat_weights$weights
  object$dat_weights$weights = NULL
  var_weights0 = var(weights0[tre_ind == 0])
  cv_weights0 = sd(weights0[tre_ind == 0])/mean(weights0[tre_ind == 0])
  eff_size0 = sum(weights0[tre_ind == 0])^2/sum(weights0[tre_ind == 0]^2)
  
  dat = object$dat_weights
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  dat[, ind] = NULL
  dat[, out] = NULL
    
  bal_cov = bal$bal_cov
  means_b0 = colMeans(dat[tre_ind == 0,])
  sds_b0 = apply(dat[tre_ind == 0,], 2, sd)
  means_b1 = colMeans(dat)
  sds_b1 = apply(dat, 2, sd)

  temp = dat
  temp[is.na(temp)] = 0
  means_a0 = as.vector(t(weights0)%*%as.matrix(temp))
  means_a1 = means_b1
  sds_a1 = sds_b1

  target = rep(NA, length(means_b0))
  tolerance = rep(NA, length(means_b0))
  
  target[match(bal_cov,colnames(dat))] = object$bal$bal_tar
  tolerance[match(bal_cov,colnames(dat))] = object$bal$bal_tol
  
  if (object$bal$bal_std == TRUE) {
    dif_b0 = as.vector(abs(target - means_b0)/sds_b0)
    dif_b1 = as.vector(abs(target - means_b1)/sds_b1)
  } else if (object$bal$bal_std == FALSE) {
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
  
  if (object$bal$bal_std == TRUE) {
    dif_a0 = as.vector(abs(target - means_a0)/sds_b0)
    dif_a1 = as.vector(abs(target - means_a1)/sds_b1)
  } else if (object$bal$bal_std == FALSE) {
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
  cat("Effective sample size of the non-missing sample: ", format(eff_size0, digits = digits),"\n")
  cat("\n")
  cat("Means of the non-missing and total sample before weighting: ", "\n")
  print(tab_b1, digits = digits)
  cat("\n")
  cat("ASDM of the non-missing and total sample before weighting: ", "\n")
  print(tab_b2, digits = digits)
  cat("\n")
  cat("Means of the non-missing and total sample after weighting: ", "\n")
  print(tab_a1, digits = digits)
  cat("\n")
  cat("ASDM of the non-missing and total sample after weighting: ", "\n")
  print(tab_a2, digits = digits)
  cat("\n")
  cat("Shadow prices for the non-missing sample: ", "\n")
  print(object$dual_table, digits = digits)
  cat("\n")

  invisible(list(variance = var_weights0, coef_variance = cv_weights0, 
                 mean_table = list(means_before = tab_b1, ASDM_before = tab_b2,
                                   means_after = tab_a1, ASDM_after = tab_a2), 
                 dual_table = object$dual_table))
}

#' Visualize output from \code{sbw}
#'
#' @description Function for summarizing the output from \code{\link[sbw]{sbw}}.
#'
#' @param object an object from the class \code{sbwcau} or \code{sbwpop} obtained after using \code{\link[sbw]{sbw}}.
#' @param digits The number of significant digits that will be displayed. The default is 6.
#' @param ... ignored arguments
#' 
#' @importFrom spatstat unnormdensity
#' 
#' @examples 
#' # Please see the examples in function sbw.
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