.sbwcaufix = function(dat, ind, out, bal, wei, sol, par) {
  if (class(dat[, ind]) == "factor") {
    dat[, ind] = as.numeric(as.character(dat[, ind]))
  } 
  if (sum(dat[, ind] != 1 & dat[, ind] != 0) > 0)
    stop(paste("Please input a binary or a logical variable for \"", ind, "\".", sep = ""))
  
  # Transform from factor to numeric
  fac_ind = sapply(dat, is.factor)
  dat[fac_ind] = lapply(dat[fac_ind], function(x) as.numeric(as.character(x)))
  # Preprocess for cate
  if(par$par_est %in% c("cate","pop")) {
    # Calculate target
    if (class(par$par_tar) == "character") {
      dat = subset(dat, eval(parse(text = par$par_tar)))
      bal$bal_tar = colMeans(as.matrix(dat[, bal$bal_cov]))
    } else if (class(par$par_tar) == "numeric") {
      if (sum(fac_ind) >= 1) {
        dat = dat[apply(dat[fac_ind] == par$par_tar[match(names(fac_ind), names(par$par_tar))][fac_ind], 1, prod, na.rm = TRUE) %in% 1,]
      }
      bal$bal_tar = par$par_tar
    } else if (class(par$par_tar) == "NULL") {
      bal$bal_tar = colMeans(as.matrix(dat[, bal$bal_cov]))
    }
  }

  if (sum(dat[ ,ind] == 0) == 0) {
    stop("Positivity is not satisfied.")
  } 
  if (sum(dat[ ,ind] == 1) == 0) {
    stop("Positivity is not satisfied.")
  } 
  
  # Order dat by ind
  ord = order(dat[, ind], decreasing = FALSE)
  dat = dat[ord, ]
  # Divide dat into list of data frames by levels of ind
  dat_level = by(dat, dat[, ind], function(x) x)
  if (par$par_est %in% c("ate", "cate")) {
    # Calculate target
    bal$bal_tar = colMeans(as.matrix(dat[, bal$bal_cov]))
    sd_target = apply(as.matrix(dat[, bal$bal_cov]), 2, sd)
    sbwfix_level = lapply(dat_level, .sbwauxfix, bal = bal, wei = wei, sol = sol, sd_target = sd_target)
    # Get weights
    weights = lapply(sbwfix_level, function(x) x$dat_weights$sbw_weights)
    # Calculate effective sample size
    effective_sample_size = lapply(weights, function(x) sum(x)^2/sum(x^2))
    # Update the outputs
    weights = unlist(weights)
    objective_value = lapply(sbwfix_level, function(x) x$objective_value)
    time = lapply(sbwfix_level, function(x) x$time)
    status = lapply(sbwfix_level, function(x) x$status)
    shadow_price = lapply(sbwfix_level, function(x) x$shadow_price)
    balance_parameters = lapply(sbwfix_level, function(x) x$balance_parameters)
  } else if (par$par_est %in% c("att", "atc", "pop")) {
    if (par$par_est %in% "att") {
      # Calculate target
      bal$bal_tar = colMeans(as.matrix(dat[which(dat[, ind] == 1), bal$bal_cov]))
      sd_target = apply(as.matrix(dat[which(dat[, ind] == 1), bal$bal_cov]), 2, sd)
      sbwfix_level = .sbwauxfix(dat_level[[1]], bal = bal, wei = wei, sol = sol, sd_target = sd_target)
      # Get weights
      weights = sbwfix_level$dat_weights$sbw_weights
      # Calculate effective sample size
      effective_sample_size = sum(weights)^2/sum(weights^2)
      # Update the data frame
      weights = c(weights, rep(1/(nrow(dat) - length(weights)), nrow(dat) - length(weights)))
    } else if (par$par_est %in% "atc") {
      bal$bal_tar = colMeans(as.matrix(dat[which(dat[, ind] == 0), bal$bal_cov]))
      sd_target = apply(as.matrix(dat[which(dat[, ind] == 0), bal$bal_cov]), 2, sd)
      sbwfix_level = .sbwauxfix(dat_level[[2]], bal = bal, wei = wei, sol = sol, sd_target = sd_target)
      # Get weights
      weights = sbwfix_level$dat_weights$sbw_weights
      # Calculate effective sample size
      effective_sample_size = sum(weights)^2/sum(weights^2)
      # Update the data frame
      weights = c(rep(1/(nrow(dat) - length(weights)), nrow(dat) - length(weights)), weights)
    } else if (par$par_est %in% "pop") {
      sd_target = apply(as.matrix(dat[, bal$bal_cov]), 2, sd)
      sbwfix_level = .sbwauxfix(dat_level[[1]], bal = bal, wei = wei, sol = sol, sd_target = sd_target)
      # Get weights
      weights = sbwfix_level$dat_weights$sbw_weights
      # Calculate effective sample size
      effective_sample_size = sum(weights)^2/sum(weights^2)
      # Update the data frame
      weights = c(weights, rep(0, nrow(dat) - length(weights)))
    }
    # Update the outputs
    objective_value = sbwfix_level$objective_value
    time = sbwfix_level$time
    status = sbwfix_level$status
    shadow_price = sbwfix_level$shadow_price
    balance_parameters = sbwfix_level$balance_parameters
  }
  # Update the data frame
  dat_weights = dat
  dat_weights$sbw_weights = weights
  dat_weights = dat_weights[order(ord), ]
  dat_weights[fac_ind] = lapply(dat_weights[fac_ind], function(x) as.factor(x))
  
  output = list(ind = ind, out = out, bal = bal, objective_value = objective_value, effective_sample_size = effective_sample_size, time = time, status = status, dat_weights = dat_weights, shadow_price = shadow_price, balance_parameters = balance_parameters, par = par)
  return(output)
}
