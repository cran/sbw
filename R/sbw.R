#' Stable balancing weights for causal contrasts and population means.
#'
#' @description Function for finding stable weights (that is, weights of minimum variance) that approximately balance the empirical distribution of the observed covariates.
#'
#' @param dat data, a data frame with a treatment assignment or missingness indicator, covariates, and possibly outcomes (which are optional).
#' @param ind treatment assignment or missingness indicator, a string with the name of the binary treatment or missingness indicator, equal to 1 if treated (missing) and 0 otherwise. 
#' When \code{par$par_est = "aux"}, \code{ind} is omitted.
#' @param out outcome, a vector of strings with the names of the outcome variables. The default is \code{NULL}.
#' @param bal balance requirements, a list with the requirements for covariate balance with the form
#' \code{list(bal_cov, bal_alg, bal_tol, bal_std, bal_gri, bal_sam)}, where:
#'
#' \code{bal_cov} balance covariates, a vector of strings with the names of the covariates in \code{dat} to be balanced.
#' In simple applications, the balance covariates in \code{bal_cov} will be the column names of \code{dat} 
#' (without including the treatment or outcome variables) for the original covariates in the data set. The covariates need to be either continuous or binary. 
#' Categorical covariates need to be transformed into dummies. In more complex applications, the covariates in \code{dat} can be 
#' transformations of the original covariates in order to balance higher order single or multidimensional moments, or other basis functions. If the transformations of the covariates are indicators of the quantiles
#' of the empirical distribution of a covariate, then balancing all these indicators will tend to balance the entire marginal distribution
#' of the covariate.
#'
#' \code{bal_alg} balance algorithm, a logical that indicates whether the tuning algorithm in Wang and Zubizarreta (2020) is 
#' to be used for automatically selecting the degree of approximate covariates balance.  The default is \code{TRUE}.  
#' See the argument \code{bal_gri} below for the candidate values for the degree of approximate covariate balance.
#' 
#' \code{bal_tol} balance tolerances, a scalar or vector of scalars
#' that define the tolerances or maximum differences in means after weighting for the covariates (or transformations thereof) defined in \code{bal_cov}.
#' Note that if \code{bal_tol} is a vector, then its length has to be equal to the length of \code{bal_cov}. 
#' Otherwise, the first element in \code{bal_tol} will be taken as the balance tolerance for all the constraints in \code{bal_cov}.
#' 
#' \code{bal_std} balance tolerances in standard deviations, a string that represent 
#' how the tolerances are adjusted. If \code{bal_std = "group"}, the tolerances 
#' are proportional to the standard deviations in the group/groups to be weighted.
#' If \code{bal_std = "target"}, the tolerances are proportional to the standard deviations 
#' in the target group. If \code{bal_std = "manual"}, the tolerances equal to \code{bal_tol}.
#' The default is \code{"group"}.
#' 
#' \code{bal_gri} grid of values for the tuning algorithm \code{bal_alg}, a vector of candidate values for the degree of approximate covariate balance. 
#' The default is \code{c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)}. 
#' The computational time is roughly proportional to the number of grid values.
#' 
#' \code{bal_sam} number of replicates to be used in \code{bal_alg}, an integer specifying the number of bootstrap sample replicates 
#' to be used to select the degree of approximate covariate balance.  See Wang and Zubizarreta (2020) for details. The default is \code{1000}. 
#' 
#' @param wei weighting constraints, a list with all the weighting constraints with the form
#' \code{list(wei_sum, wei_pos)}, where:
#' 
#' \code{wei_sum} sum of weights, a logical variable indicating whether the weights are constrained to sum up to one, or whether their sum 
#' is unconstrained. The default is \code{TRUE} for the sum of weights equal to one. Note that if \code{wei_sum = TRUE}, then \code{wei_pos = TRUE}.
#' 
#' \code{wei_pos} positive or zero (non-negative) weights, a logical variable indicating whether the weights are constrained to be non-negative, or whether they
#' are unconstrained. The default is \code{TRUE} for non-negative weights.  Again, note that if \code{wei_sum = TRUE}, then \code{wei_pos = TRUE}.
#'
#' @param sol solver, a list that specifies the solver option with the form
#' 
#' \code{list(sol_nam, sol_dis, sol_pog)}, where:
#' 
#' \code{sol_nam} solver name, a string equal to either \code{"cplex"}, \code{"gurobi"}, \code{"mosek"}, \code{"osqp"}, \code{"pogs"}, or \code{"quadprog"}.
#' CPLEX, \href{https://www.gurobi.com/documentation/current/refman/r_ins_the_r_package.html}{Gurobi} and MOSEK are commercial solvers, but free for academic users. 
#' \href{http://foges.github.io/pogs/stp/r}{POGS} and QUADPROG are free for all. In our experience, POGS is the fastest solver option
#' and able to handle larger datasets, but it can be difficult to install for non-Mac users 
#' and more difficult to calibrate. MOSEK is more stable than POGS and faster. 
#' The default option is \code{sol_nam = "quadprog"}. 
#'
#' \code{sol_dis} solver display, a logical variable indicating whether the output is to be displayed or not.
#' The default is \code{FALSE}. This option is specific to \code{"cplex"}, \code{"gurobi"}, \code{"mosek"}, \code{"pogs"}, and \code{"osqp"}.
#'
#' \code{sol_pog} solver options specific to \code{"pogs"}, with the following default parameters:
#'
#' \code{sol_pog = list(sol_pog_max_iter = 100000, sol_pog_rel_tol = 1e-4,} 
#' 
#' \code{sol_pog_abs_tol = 1e-4, sol_pog_gap_stp = TRUE, sol_pog_adp_rho = TRUE)}.
#' 
#' See the POGS manual for details.
#'  
#' @param par parameter of interest, a list describing the parameter of interest or estimand with the form
#' \code{list(par_est, par_tar)}, where:
#' 
#' \code{par_est} estimand. For causal inference, a string equal to:
#' \code{"att"} (Average Treatment effect on the Treated), 
#' \code{"atc"} (Average Treatment effect on the Controls), 
#' \code{"ate"} (Average Treatment Effect), 
#' \code{"cate"} (Conditional Average Treatment Effect). 
#' For estimation with incomplete outcome data, a string equal to:
#' \code{"pop"} (General population means) or
#' \code{"aux"} (Means for a population specified by the user). The default is \code{"att"}.
#' 
#' \code{par_tar} target, a string, or a vector of scalars. 
#' It specifies the targeted population for inference in terms of the observed covariates 
#' when \code{par_est = "cate"}, \code{"pop"} or \code{"aux"}. Please see the examples. 
#' 
#' @param mes a logical variable indicating whether the messages are printed.
#' 
#' @import quadprog slam
#' @importFrom MASS mvrnorm
#' @importFrom Matrix sparseMatrix Diagonal bdiag t
#' 
#' @return A list with the following elements:
#' @return \code{dat_weights}, a data frame with the optimal weights \code{dat_weights$sbw_weights};
#' @return \code{ind}, an argument provided by the user;
#' @return \code{out}, an argument provided by the user;
#' @return \code{bal}, an argument provided by the user;
#' @return \code{wei}, an argument provided by the user;
#' @return \code{sol}, an argument provided by the user;
#' @return \code{par}, an argument provided by the user;
#' @return \code{effective_sample_size}, effective sample size/sizes for the weighted group/groups;
#' @return \code{objective_value}, value/values of the objective function/functions at the optimum;
#' @return \code{status}, status of the solution. If the optimal weights are found, \code{status = optimal};
#' otherwise, the solution may be not optimal or not exist, in which case an error will be returned with details specific to the solver used.
#' For the solver "quadprog", the status code is missing, therefore, \code{status = NA};
#' @return \code{time}, time elapsed to find the optimal solution;
#' @return \code{shadow_price}, dual variables or shadow prices of the covariate balance constraints;
#' @return \code{balance_parameters}, details of the balance parameters;
#' @return \code{cstat}, covariate balance statistic used in Wang and Zubizarreta (2020).
#' A magnitude to be minimized to select the degree of approximate balance in \code{bal$bal_gri}.
#' 
#' @source{\url{https://www.ibm.com/products/ilog-cplex-optimization-studio}} 
#' @source{\url{https://www.gurobi.com/products/gurobi-optimizer/}}
#' @source{\url{https://www.mosek.com/products/mosek/}} 
#' @source{\url{http://foges.github.io/pogs/stp/r}} 
#' 
#' @references Chattopadhyay, A., Hase, C. H., and Zubizarreta, J. R. (2020), "Balancing Versus Modeling Approaches to Weighting in Practice," \emph{Statistics in Medicine}, 39, 3227-3254.
#' @references Kang, J. D. Y., and Schafer, J. L. (2007), "Demystifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data," \emph{Statistical Science}, 22, 523-539.
#' @references Stuart, E. A. Matching methods for causal inference: a review and a look forward. Statistical Science 2010; 25(1): 1-21.
#' @references Wang, Y., and Zubizarreta, J. R. (2020), "Minimal Dispersion Approximately Balancing Weights: Asymptotic Properties and Practical Considerations," \emph{Biometrika}, 107, 93-105.
#' @references Zubizarreta, J. R. (2015), "Stable Weights that Balance Covariates for Estimation with Incomplete Outcome Data," \emph{Journal of the American Statistical Association}, 110, 910-922.
#' 
#' @examples
#' # Simulate data
#' kangschafer = function(n_obs) {
#'  # Z are the true covariates
#'  # t is the indicator for the respondents (treated)
#'  # y is the outcome
#'  # X are the observed covariates
#'  # Returns Z, t y and X sorted in decreasing order by t
#'  Z = MASS::mvrnorm(n_obs, mu=rep(0, 4), Sigma=diag(4))
#'  p = 1/(1+exp(Z[, 1]-.5*Z[, 2]+.25*Z[, 3]+.1*Z[, 4]))
#'  t = rbinom(n_obs, 1, p)
#'  Zt = cbind(Z, p, t)
#'  Zt = Zt[order(t), ]
#'  Z = Zt[, 1:4]
#'  p = Zt[, 5]
#'  t = Zt[, 6]
#'  y = 210+27.4*Z[, 1]+13.7*Z[, 2]+13.7*Z[, 3]+13.7*Z[, 4]+rnorm(n_obs)
#'  X = cbind(exp(Z[, 1]/2), (Z[, 2]/(1+exp(Z[, 1])))+10, (Z[, 1]*Z[, 3]/
#' 25+.6)^3, (Z[, 2]+Z[, 4]+20)^2)
#'  return(list(Z=Z, p=p, t=t, y=y, X=X))
#' }
#' set.seed(1234)
#' n_obs = 200
#' aux = kangschafer(n_obs)
#' Z = aux$Z
#' p = aux$p
#' t = aux$t
#' y = aux$y
#' X = aux$X
#'
#' # Generate data frame
#' t_ind = t
#' bal_cov = X
#' data_frame = as.data.frame(cbind(t_ind, bal_cov, y))
#' names(data_frame) = c("t_ind", "X1", "X2", "X3", "X4", "Y")
#'
#' # Define treatment indicator and 
#' t_ind = "t_ind"
#' # moment covariates
#' bal = list()
#' bal$bal_cov = c("X1", "X2", "X3", "X4")
#'
#' # Set tolerances
#' bal$bal_tol = 0.02
#' bal$bal_std = "group"
#' 
#' # Solve for the Average Treatment Effect on the Treated, ATT (default)
#' bal$bal_alg = FALSE
#' sbwatt_object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal)
#' 
#' # # Solve for a Conditional Average Treatment Effect, CATE
#' # sbwcate_object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "cate", par_tar = "X1 > 1 & X3 <= 0.22"))
#' 
#' # # Solve for the population mean, POP
#' # tar = colMeans(bal_cov)
#' # names(tar) = bal$bal_cov
#' # sbwpop_object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "pop"))
#' 
#' # # Solve for a target population mean, AUX
#' # sbwaux_object = sbw(dat = data_frame, bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "aux", par_tar = tar*1.05))
#' 
#' # # Solve for the ATT using the tuning algorithm
#' # bal$bal_alg = TRUE
#' # bal$bal_sam = 1000
#' # sbwatttun_object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "att", par_tar = NULL))
#' 
#' # Check
#' summarize(sbwatt_object)
#' # summarize(sbwcate_object)
#' # summarize(sbwpop_object)
#' # summarize(sbwaux_object)
#' # summarize(sbwatttun_object)
#' 
#' # Estimate
#' estimate(sbwatt_object)
#' # estimate(sbwcate_object)
#' # estimate(sbwpop_object)
#' # estimate(sbwatttun_object)
#' 
#' # Visualize
#' visualize(sbwatt_object)
#' # visualize(sbwcate_object)
#' # visualize(sbwpop_object)
#' # visualize(sbwaux_object)
#' # visualize(sbwatttun_object)
#' 
#' @export
#'
sbw = function(dat, ind = NULL, out = NULL, bal = list(bal_cov, bal_alg = TRUE, bal_tol, bal_std = "group", bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1), bal_sam = 1000), wei = list(wei_sum = TRUE, wei_pos = TRUE), sol = list(sol_nam = "quadprog", sol_dis = FALSE), par = list(par_est = "att", par_tar = NULL), mes = TRUE) {
  if (is.null(bal$bal_alg)) {
    bal$bal_alg = TRUE
  } else if (!is.logical(bal$bal_alg)) {
    stop("Please assign a logical to bal$bal_alg.")
  } 
  
  if (is.null(bal$bal_std)) {
    bal$bal_std = "group"
  } 
  if (!(bal$bal_std %in% c("group", "target", "manual"))) {
    stop("bal$bal_std should be equal to one of 'group', 'target', 'manual'.")
  }
  
  if (is.null(bal$bal_gri)) {
    bal$bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)
  } else if (!is.numeric(bal$bal_gri)) {
    stop("Please assign a numeric vector to bal$bal_gri.")
  } 

  if (is.null(bal$bal_sam)) {
    bal$bal_sam = 1000
  } else if (!is.numeric(bal$bal_sam)) {
    stop("Please assign an integer to bal$bal_sam.")
  }
  
  if (!par$par_est %in% c("aux")) {
    if (sum(1 - is.na(match(ind, colnames(dat)))) == 0) {
      stop("Please specify a correct string for ind.")
    }
  }
  
  if ("sbw_weights" %in% colnames(dat)) {
    warning("The 'sbw_weights' column in 'dat' will be replaced by the optimal weights. Please change the name of the original 'sbw_weights' to avoid unexpected replacement.")
  }
  
  if (!par$par_est %in% c("aux")) {
    if (bal$bal_alg == FALSE) {
      if (mes == FALSE) {
        invisible(capture.output(object <- .sbwfix(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)))
      } else object = .sbwfix(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)
    } else if (bal$bal_alg == TRUE) {
      bal$bal_tol = bal$bal_gri
      if (mes == FALSE) {
        invisible(capture.output(object <- .sbwtun(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)))
      } else object = .sbwtun(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)
    }
  } else if (par$par_est %in% c("aux")) {
    if (bal$bal_std %in% "target") {
      warning("The tolerances will be 0.")
    }
    sd_target = 0
    if (bal$bal_alg == FALSE) {
      bal$bal_tar = par$par_tar
      if (mes == FALSE) {
        invisible(capture.output(object <- .sbwauxfix(dat = dat, bal = bal, wei = wei, sol = sol, sd_target = sd_target)))
      } else object = .sbwauxfix(dat = dat, bal = bal, wei = wei, sol = sol, sd_target = sd_target)
    } else if (bal$bal_alg == TRUE) {
      bal$bal_tar = par$par_tar
      if (mes == FALSE) {
        invisible(capture.output(object <- .sbwauxtun(dat = dat, bal = bal, wei = wei, sol = sol, sd_target = sd_target)))
      } else object = .sbwauxtun(dat = dat, bal = bal, wei = wei, sol = sol, sd_target = sd_target)
    }
    object$ind = ind
    object$out = out
    object$par = par
  }
  object_class = class(object)
  object_name = factor(names(object),
                       levels = c("dat_weights",
                                  "ind", "out", "bal", "wei", "sol", "par",
                                  "effective_sample_size", "objective_value", "status", "time",
                                  "shadow_price", "balance_parameters", "cstat"))
  object = object[order(object_name)]
  class(object) = object_class
  return(object)
}
