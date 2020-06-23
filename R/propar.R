if(getRversion() >= "2.15.1")
  utils::globalVariables(
    c("abline", "boxplot", "bal", "bal_cov", "bal_tar0", "bal_tar1", "bal_tol", "cor", "dat_weights", "density",
      "digits", "shadow_price", "shadow_price0", "shadow_price1", "effective_sample_size", "gri",
      "lb", "legend", "lines", "objective_value", "par", "plot", "pt", 
      "Qmat", "qt", "sam", "sd", "segments", "sense", "status",
      "tail", "target", "target0", "target1", "time", "ub",
      "var_type", "var", "weights", "weighted.mean", "..."))

# Translate the balancing problem to an optimization problem and acquire parameters for solvers.
.problemparameters = function(dat, nor, bal, normalize, w_min, sd_target) {
  # Check bal arguments
  if (length(bal$bal_cov) == 0) {
    stop("bal_cov should not be empty.")
  }
  if (length(bal$bal_tar) == 0) {
    stop("bal_tar should not be empty.")
  }
  if (length(bal$bal_tar) != length(bal$bal_cov)) {
    stop("bal_cov and bal_tar should have equal length.")
  } 
  if (length(bal$bal_tol) == 0) {
    stop("bal_tol should not be empty.")  
  }
  
  # Transform the inputs with factor class
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  
  # Get numeric bal arguments
  bal_cov = dat[, bal$bal_cov]
  bal_cov = as.matrix(bal_cov)
  bal_tar = bal$bal_tar
  names(bal_tar) = bal$bal_cov
  bal_tol = bal$bal_tol
  bal_std = ifelse(is.null(bal$bal_std), "group", bal$bal_std)
  
  if (length(bal_tol) > 1 & (length(bal$bal_cov) != length(bal_tol))) {
    warning("The length of bal$bal_cov and bal$bal_tol are different. The first element of bal$bal_tol will be used as bal$bal_tol.")
    bal_tol = rep(bal_tol[1], length(bal$bal_cov))
    names(bal_tol)  = bal$bal_cov
  }
  if (length(bal_tol) == 1) {
    bal_tol = rep(bal_tol, length(bal$bal_cov))
    names(bal_tol)  = bal$bal_cov
  }
  
  if (bal_std %in% "group") {
    bal_tol = apply(bal_cov, 2, sd)*bal_tol
  } else if (bal_std %in% "target") {
    bal_tol = sd_target*bal_tol
  } else if (!(bal_std %in% c("group", "target", "manual"))) {
    stop("bal$bal_std should be equal to one of 'group', 'target', 'manual'.")
  }
  
  # Get bal_new
  bal_new = list(bal_tar = bal_tar, bal_tol = bal_tol, bal_std = bal_std)
  
  # Get n, p
  n = nrow(bal_cov)
  p = ncol(bal_cov)
  
  # Get nor
  if (!(nor %in% c("l_2", "l_1", "l_inf"))) stop("Please use one of l_2(default), l_1 or l_inf norms.")
  
  # Build parameter cvec
  Qmat = NULL
  if (nor == "l_1") {
    cvec = c(rep(0, n), rep(1, n))
  }
  if (nor == "l_inf") {
    cvec = c(rep(0, n), 1)
  }
  if (nor == "l_2") {
    cvec = rep(0, n)
  }
  
  # Build parameter Amat
  row_ind_cur = 0
  if (nor == "l_1") {
    row_ind_nor = sort(rep(1:(n*2), n + 1)) + row_ind_cur
    col_ind_nor = matrix(rep(1:n, n*2), nrow = n)
    col_ind_nor = rbind(col_ind_nor, sort(rep((1:n) + n, 2)))
    col_ind_nor = as.vector(col_ind_nor)
    vals_l_1 = rep(c(rep(-1/n, n), -1, rep(1/n, n), -1), n)
    aux_ind = ((sort(rep(1:n, 2))-1)*(n + 1)) + sort(rep((0:(n - 1)), 2)) + c(1, 1, rep(0, (n*2)-2))
    vals_l_1[aux_ind] = rep(c(1, -1), n) + vals_l_1[aux_ind]
    row_ind_cur = max(row_ind_nor)
  }
  if (nor == "l_inf") {
    row_ind_nor = sort(rep(1:(n*2), n + 1)) + row_ind_cur
    col_ind_nor = matrix(rep(1:n, n*2), nrow = n)
    col_ind_nor = rbind(col_ind_nor, rep(n + 1, n*2))
    col_ind_nor = as.vector(col_ind_nor)
    vals_l_1 = rep(c(rep(-1/n, n), -1, rep(1/n, n), -1), n)
    aux_ind = ((sort(rep(1:n, 2))-1)*(n + 1)) + sort(rep((0:(n - 1)), 2)) + c(1, 1, rep(0, (n*2)-2))
    vals_l_1[aux_ind] = rep(c(1, -1), n) + vals_l_1[aux_ind]
    row_ind_cur = max(row_ind_nor)
  }
  row_ind_mom = sort(rep(1:(p*2), n)) + row_ind_cur
  col_ind_mom = rep(1:n, p*2)
  vals_mom = cbind(bal_cov, -bal_cov)[, sort(rep(1:p, 2)) + rep(c(0, p), p)]
  row_ind_cur = max(row_ind_mom) + 1
  
  # Add additional constraints
  if (normalize == 1) {
    row_ind_normalize = rep(row_ind_cur, n)
    col_ind_normalize = 1:n
    vals_normalize = rep(1, n)
    row_ind_cur = max(row_ind_normalize) + 1
  }
  if (nor == "l_1" | nor == "l_inf") {
    row_ind = c(row_ind_nor, row_ind_mom)
    col_ind = c(col_ind_nor, col_ind_mom)
    vals = c(vals_l_1, vals_mom)
    if (normalize == 1) {
      row_ind = c(row_ind_nor, row_ind_mom, row_ind_normalize)
      col_ind = c(col_ind_nor, col_ind_mom, col_ind_normalize)
      vals = c(vals_l_1, vals_mom, vals_normalize)
    }
  }
  if (nor == "l_2") {
    row_ind = c(row_ind_mom)
    col_ind = c(col_ind_mom)
    vals = c(vals_mom)
    if (normalize == 1) {
      row_ind = c(row_ind_mom, row_ind_normalize)
      col_ind = c(col_ind_mom, col_ind_normalize)
      vals = c(vals_mom, vals_normalize)
    }
  }
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  Amat = slam::simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  # Build parameter bvec and sort with the order of covariates and upper/lower bound
  if (nor == "l_1" | nor == "l_inf") {
    bvec = rep(0, n*2)
    bvec_mom = c(bal_tar + bal_tol, - bal_tar + bal_tol)[sort(rep(1:p, 2)) + rep(c(0, p), p)]
    bvec = c(bvec, bvec_mom)
    if (normalize == 1) {
      bvec = c(bvec, 1)
    }
  }
  
  if (nor == "l_2") {
    bvec = c(bal_tar + bal_tol, -bal_tar + bal_tol)[sort(rep(1:p, 2)) + rep(c(0, p), p)]
    if (normalize == 1) {
      bvec = c(bvec, 1)
    }
  }
  
  # Build parameters lower bound and upper bound
  if (nor == "l_1") {
    lb = rep(w_min, n)
    lb = c(lb, rep(0, n))
  }
  if (nor == "l_inf") {
    lb = rep(w_min, n)
    lb = c(lb, 0)
  }
  if (nor == "l_2") {
    lb = rep(w_min, n)
  }
  if (nor == "l_1") {
    ub = rep(Inf, n*2)
  }
  if (nor == "l_inf") {
    ub = rep(Inf, 1 + n)
  }
  if (nor == "l_2") {
    ub = rep(Inf, n)
  }
  
  # Build parameters sense
  if (nor == "l_1") {
    sense = rep("L", n*2)
    sense = c(sense, rep("L", p*2))
    if (normalize == 1) {
      sense = c(sense, "E")
    }
  }
  if (nor == "l_inf") {
    sense = rep("L", n*2)
    sense = c(sense, rep("L", p*2))
    if (normalize == 1) {
      sense = c(sense, "E")
    }
  }
  if (nor == "l_2") {
    sense = rep("L", p*2)
    if (normalize == 1) {
      sense = c(rep("L", p*2), "E")
    }
  }
  
  # Build parameters var_type
  if (nor == "l_1") {
    var_type = rep("C", n*2)
  }
  if (nor == "l_inf") {
    var_type = rep("C", 1+n)
  }
  if (nor == "l_2") {
    var_type = rep("C", n)
  }
  
  return(list(normalize = normalize, nor = nor, n = n, Amat = Amat, bvec = bvec, cvec = cvec, lb = lb, ub = ub, sense = sense, var_type = var_type, bal_new = bal_new))
}