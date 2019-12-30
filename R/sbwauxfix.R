.sbwauxfix = function(dat, bal, wei, sol, ...) {
  # Transform some arguments to inner functions
  if (wei$wei_sum == TRUE & wei$wei_pos == TRUE) {
    normalize = 1
    w_min = 0
  } else if (wei$wei_sum == FALSE & wei$wei_pos == TRUE) {
    normalize = 0
    w_min = 0
  } else if (wei$wei_sum == FALSE & wei$wei_pos == FALSE) {
    normalize = 0
    w_min = -Inf
  } else {stop("If sum is TRUE, then pos is supposed to be TRUE.")}
  nor = "l_2"
  
  # Check errors:
  # unmatched covariates
  check_cov = bal$bal_cov[is.na(match(bal$bal_cov, colnames(dat)))]
  if (length(check_cov) > 0) stop(paste(paste(check_cov, collapse = ", "), "are not found in the dat."))
  # missing values
  if (sum(is.na(dat[, bal$bal_cov])) > 0) {
    mis_value = colSums(is.na(dat[, bal$bal_cov]))
    stop(paste(paste(names(which(mis_value != 0)), collapse = ", "), "have missing values."))
  }
  # about bal, bal$bal_tar is the target
  if (length(bal$bal_cov) != length(bal$bal_tar)) stop("bal$bal_cov and par$par_tar should have the same length as well as the same order.")
  if (!is.numeric(bal$bal_tol)) {
    stop("bal$bal_tol should be numeric.")
  } else if (sum(bal$bal_tol <= 0) > 0) stop("bal$bal_tol should be positive.")
  
  # Set up optimization question
  problemparameters.object = .problemparameters(dat = dat, nor = nor, bal = bal, normalize = normalize, w_min = w_min)
  
  # Robust to capital input for sol_nam
  sol$sol_nam = tolower(sol$sol_nam)
  
  # Choose solvers
  if (sol$sol_nam == "cplex") {
    if (requireNamespace('Rcplex', quietly = TRUE)) {
      trace = ifelse(is.null(sol$sol_dis), 0, as.numeric(sol$sol_dis))
      ptm = proc.time()
      sbw.object = .sbwpricplex(problemparameters.object, trace = trace) 
      time = (proc.time()-ptm)[3]
    }
  }
  else if (sol$sol_nam == "pogs") {
    if (requireNamespace("pogs", quietly=TRUE)) {
      # Set defaults
      max_iter = ifelse(is.null(sol$sol_pog$sol_pog_max_iter), 100000, sol$sol_pog$sol_pog_max_iter)
      rel_tol = ifelse(is.null(sol$sol_pog$rel_tol), 1e-04, sol$sol_pog$rel_tol)
      abs_tol = ifelse(is.null(sol$sol_pog$abs_tol), 1e-04, sol$sol_pog$abs_tol)
      gap_stop = ifelse(is.null(sol$sol_pog$gap_stop), TRUE, sol$sol_pog$gap_stop)
      adaptive_rho = ifelse(is.null(sol$sol_pog$adaptive_rho), TRUE, sol$sol_pog$adaptive_rho)
      verbose = ifelse(is.null(sol$sol_dis), 0, as.numeric(sol$sol_dis))
      params = list(rel_tol = rel_tol, abs_tol = abs_tol, max_iter = max_iter, adaptive_rho = adaptive_rho, verbose = verbose, gap_stop = gap_stop)
      ptm = proc.time()
      sbw.object = .sbwpripogs(problemparameters.object, params) 
      time = (proc.time()-ptm)[3]
    }
  }
  else if (sol$sol_nam == "quadprog") {
    ptm = proc.time()
    sbw.object = .sbwpriquadprog(problemparameters.object) 
    time = (proc.time()-ptm)[3]
  }
  else if (sol$sol_nam == "gurobi") {
    if (requireNamespace('gurobi', quietly=TRUE)) {
      sol_dis = ifelse(is.null(sol$sol_dis), 0, as.numeric(sol$sol_dis))
      params = list(OutputFlag = sol_dis)
      ptm = proc.time()
      sbw.object = .sbwprigurobi(problemparameters.object, params = params) 
      time = (proc.time()-ptm)[3]
    }
  }
  else if (sol$sol_nam == "mosek") {
    if (requireNamespace('Rmosek', quietly=TRUE)) {
      verbose = ifelse(is.null(sol$sol_dis), 0, as.numeric(paste(as.numeric(sol$sol_dis), 0, sep = "")))
      ptm = proc.time()
      sbw.object = .sbwprimosek(problemparameters.object, verbose = verbose) 
      time = (proc.time()-ptm)[3]
    }
  } else {stop("The assignment of solver is not found, please choose one of the available solvers.")}

  # Save target for evaluation functions 
  target = problemparameters.object$bal_new
  if (length(bal$bal_cov) != length(bal$bal_tol)) {
    bal$bal_tol = bal$bal_tol[1]
  }
  target$bal_tol_ori = bal$bal_tol
  rm(problemparameters.object)
  
  dat_weights = dat
  dat_weights$weights = as.numeric(as.character(sbw.object$weights))
  eff_size = sum(dat_weights$weights)^2/sum(dat_weights$weights^2)
  
  output = list(bal = bal, wei = wei, sol = sol, obj_total = sbw.object$obj_total, eff_size = eff_size, time = time, status = sbw.object$status, dat_weights = dat_weights, dual_table = sbw.object$dual_table, target = target)
  class(output) = "sbwaux"
  return(output)
}