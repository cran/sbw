.sbwauxtun = function(dat, bal, wei, sol, run = TRUE, ...) {
  # Check errors:
  # about bal, bal$bal_tar is the target
  if (length(bal$bal_cov) != length(bal$bal_tar)) stop("bal$bal_cov and par$par_tar should have the same length as well as the same order.")
  if (!is.numeric(bal$bal_gri)) {
    stop("bal$bal_gri should be a numeric vector.")
  } else if (sum(bal$bal_gri <= 0) > 0) stop("bal$bal_gri should be all positive values.")
  
  Cstat = function(tol) {
    bal$bal_tol = tol
    outSBWAUX = .sbwauxfix(dat = dat, bal = bal, wei = wei, sol = sol)
    weights = outSBWAUX$dat_weights$weights
    n = length(weights)
    C_k = 0
    for (k in 1:bal$bal_sam) {
      ind = sample(n, n, replace = TRUE)
      if (bal$bal_std == TRUE) {
        C_k = C_k + mean(abs((t(weights[ind])%*%as.matrix(dat[ind, bal$bal_cov])/sum(weights[ind]) - bal$bal_tar)/apply(as.matrix(dat[, bal$bal_cov]), 2, sd)))
      } else if (bal$bal_std == FALSE) {
        C_k = C_k + mean(abs((t(weights[ind])%*%as.matrix(dat[ind, bal$bal_cov])/sum(weights[ind]) - bal$bal_tar)))
      }
    }
    C_k/bal$bal_sam
  }
  Cstat.object = lapply(bal$bal_gri, Cstat)
  bal$bal_tol = bal$bal_gri[which.min(Cstat.object)]
  cstat = unlist(Cstat.object)
  if (run == TRUE) {
    sbwfix.object = .sbwauxfix(dat = dat, bal = bal, wei = wei, sol = sol)
    output = list(cstat = cstat, bal = bal, wei = wei, sol = sol, 
                  obj_total = sbwfix.object$obj_total, eff_size = sbwfix.object$eff_size, 
                  time = sbwfix.object$time, status = sbwfix.object$status, 
                  dat_weights = sbwfix.object$dat_weights, 
                  dual_table = sbwfix.object$dual_table, target = sbwfix.object$target)
  } else {
    output = list(cstat =  cstat)
  }
  class(output) = "sbwaux"
  return(output)
}
