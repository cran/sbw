.sbwtun = function(dat, ind, out, bal, wei, sol, par, ...) {
  # Check errors
  check_cov = bal$bal_cov[is.na(match(bal$bal_cov, colnames(dat)))]
  if (length(check_cov) > 0) stop(paste(paste(check_cov, collapse = ", "), "are not found in the dat."))
  if (sum(is.na(dat[, bal$bal_cov])) > 0) {
    mis_value = colSums(is.na(dat[, bal$bal_cov]))
    stop(paste(paste(names(which(mis_value != 0)), collapse = ", "), "have missing values."))
  }
  if (class(par) == "list") {
    est = par$par_est
  } else {stop("Please input a list for argument \"par\".")}
  
  if (est %in% c("att", "atc", "ate", "cate", "pop")) {
    output = .sbwcautun(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)
  } else {stop("Please input one of \"att\", \"atc\", \"ate\", \"cate\", \"pop\" for argument est.")}
  
  bal$bal_tol = output$bal$bal_tol
  bal$bal_tar = output$bal$bal_tar
  
  output = list(ind = ind, out = out,  bal = bal, cstat = output$cstat, wei = wei, sol = sol, par = par, objective_value = output$objective_value, effective_sample_size = output$effective_sample_size, time = output$time, status = output$status, dat_weights = output$dat_weights, shadow_price = output$shadow_price, balance_parameters = output$balance_parameters)
  
  if (est %in% c("att", "atc", "ate", "cate")) {
    class(output) = "sbwcau"
  } else if (est %in% c("pop")) {
    class(output) = "sbwpop"
  }
  output
}