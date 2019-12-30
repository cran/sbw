#' Stable balancing weights for causal contrasts and population means.
#'
#' @description Function for finding stable weights (this is, weights of minimum variance) that approximately balance the empirical distribution of the observed covariates.
#'
#' @param dat data, a data frame including a treatment or missingness indicator plus covariates; outcomes are optional.
#' @param ind treatment or missingness indicator, a string with the name of the binary treatment or missingness, equal to 1 if treated (missing) and 0 otherwise. 
#' When \code{par$par_est = "aux"}, \code{ind} is omitted.
#' @param out outcome, a vector of strings with the names of the outcome variables. The default is \code{NULL}.
#' @param bal balance requirements, a list with requirements for covariate balance with the form
#' \code{list(bal_cov, bal_alg, bal_tol, bal_std, bal_gri, bal_sam)}, where:
#'
#' \code{bal_cov}: balance covariates, a vector of strings with the names of the covariates in \code{dat} to be balanced.
#' In simple applications, the balance covariates in \code{bal_cov} will be the column names of \code{dat} (of course,
#' without including the treatment or outcome variables) for the original covariates in the data set. The covariates need to be either continuous or binary. 
#' Categorical variables need to be transformed into dummy variables. In more complex applications, the covariates in \code{dat} can be 
#' transformations of the original covariates in order to balance higher order single dimensional moments such as variances and skewness,
#' and multidimensional moments such as correlations. If the transformations of the covariates are indicators of the quantiles
#' of the empirical distribution of a covariate, then balancing all these indicators will tend to balance the entire marginal distribution
#' of the covariate.
#'
#' \code{bal_alg} balance algorithm, a logical that indicates whtether the tuning algorithm in Wang and Zubizarreta (2019) is 
#' to be used for automatically selecting the degree of approximate covariates balance.  The default is TRUE.  
#' See the argument \code{bal_gri} below for the candidate values for the degree of approximate covariate balance.
#' 
#' \code{bal_tol}: balance tolerances, a scalar or vector of scalars
#' defining the tolerances or maximum differences in means after weighting for the covariates defined in \code{bal_cov}.
#' Note that if \code{bal_tol} is a vector then its length has to be equal to the length of \code{bal_cov}. 
#' Otherwise, the first element in the vector will be taken as the balance tolerance for all the constraints in \code{bal_cov}.
#' 
#' \code{bal_std}: balance tolerances in standard deviations, a logical that indicates
#' whether the tolerances specified in \code{bal_tol} are expressed in the original units of the covariates
#' or in standard deviations. The default is \code{TRUE}, meaning that the tolerances are expressed in standard deviations.
#' 
#' \code{bal_gri}: grid of values for the tuning algorithm \code{bal_alg}, a vector of candidate values for the degree of approximate covariate balance. 
#' The default is \code{c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1)}. 
#' The computational time is roughlly proportional to the number of grid values.
#' 
#' \code{bal_sam}: number of replicates to be used in \code{bal_alg}, an integer specifying the number of bootstrap sample replicates 
#' to be used to select the degree of approximate covariate balance.  See Wang and Zubizarreta (2019) for details. The default is \code{1000}. 
#' 
#' @param wei weighting constraints, a list with all the weighting constraints with the form
#' \code{list(wei_sum, wei_pos)}, where:
#' 
#' \code{wei_sum}: sum of weights, a logical variable indicating whether the weights are constrained to sum up to one, or their sum 
#' is unconstrained. The default is \code{TRUE} for the sum of weights equal to one. Note that if \code{wei_sum = TRUE}, then \code{wei_pos = TRUE}.
#' 
#' \code{wei_pos}: positive or zero (non-negative) weights, a logical variable indicating whether the weights are constrained to be non-negative, or they
#' are unconstrained. The default is \code{TRUE} for non-negative weights.  Again, note that if \code{wei_sum = TRUE}, then \code{wei_pos = TRUE}.
#'
#' @param sol solver, a list specifying the solver option with the form
#' \code{list(sol_nam, sol_dis, sol_pog)} where:
#' 
#' \code{sol_nam}: solver name, a string equals to one of "cplex", "gurobi", "mosek", "pogs", "quadprog".
#' CPLEX, Gurobi and MOSEK are commercial solvers, but free for academic users. 
#' POGS and QUADPROG are free for all. In our experience, POGS is the fastest solver option
#' and able to handle larger datasets, but it can be difficult to install for non-Mac users 
#' and more difficult to callibrate. MOSEK is more stable than POGS and fast. 
#' The default option is solver = "quadprog". 
#'
#' \code{sol_dis}: solver display, a logical variable indicating whether the output is to be displayed or not.
#' The default is FALSE. This option is specific to "cplex", "gurobi", "mosek" and "pogs".
#'
#' \code{sol_pog}: solver options specific to "pogs", with the following default parameters:
#'
#' \code{sol_pog = list(sol_pog_max_iter = 100000, sol_pog_rel_tol = 1e-4,} 
#' \code{sol_pog_abs_tol = 1e-4, sol_pog_gap_stp = TRUE, sol_pog_adp_rho = TRUE)}.
#' 
#' See the POGS manual for details.
#' @source{\url{http://foges.github.io/pogs/stp/r}} 
#'  
#' @param par parameter of interest, a list describing the parameter of interest or estimand with the form
#' \code{list(par_est, par_tar)}, where
#' 
#' \code{par_est} estimand, a string. For causal inference, a string equals to one of:
#' "att" (Average Treatment effect among the Treated), 
#' "atc" (Average Treatment effect among the Controls), 
#' "ate" (Average Treatment Effect), 
#' "cate" (Conditional Average Treatment Effect). 
#' For estimation with incomplete outcome data, a string equal to:
#' "pop" (Population Means) or
#' "aux" (Specified Targets by users). The default is "att".
#' 
#' \code{par_tar} target, a string or a vector that specifies the targeted population for inference in terms of 
#' the observed covariates if \code{par_est} = "cate", "pop" or "aux" with the form 
#' of determine statements. Please see the examples. It also accepts a numeric vector 
#' as specified balance targets for \code{bal$bal_cov}.
#' 
#' @import quadprog slam
#' @importFrom MASS mvrnorm
#' @importFrom Matrix sparseMatrix Diagonal bdiag t
#' 
#' @return A list with the following elements:
#' @return \code{ind}{, an inputted argument}
#' @return \code{out}{, an inputted argument}
#' @return \code{bal}{, an inputted argument}
#' @return \code{wei}{, an inputted argument}
#' @return \code{sol}{, an inputted argument}
#' @return \code{par}{, an inputted argument}
#' @return \code{obj_total}{, value/values of the objective function/functions at the optimum;}
#' @return \code{eff_size}{, effective sample size/sizes for the weighted group/groups;}
#' @return \code{time}{, time elapsed to find the optimal solution;}
#' @return \code{status}{, status of the solution. If the optimal weights are found,} \code{status = optimal}
#' {; otherwise, the solution may be not optimal or not exist, in which case an error will be returned with details specific to the solver used.}
#' { For the solver "quadprog", the status code is missing, therefore, } \code{status = NA} {;}
#' @return \code{dat_weights}{, data frame with the optimal weights named by weights;}
#' @return \code{dual_table}{, dual variables or shadow prices of the covariate balancing constraints;}
#' @return \code{target}{, details of the balance targets, which are saved for evaluation uses.}
#' 
#' @references Chattopadhyay, A., Hase, C. H., and Zubizarreta, J. R. (2019), "Balancing Versus Modeling Approaches to Weighting in Practice," submitted.
#' @references Kang, J. D. Y., and Schafer, J. L. (2007), "Demistifying Double Robustness: A Comparison of Alternative Strategies for Estimating a Population Mean from Incomplete Data," Statistical Science, 22, 523-539.
#' @references Wang, Y., and Zubizarreta, J. R. (2019), "Minimal Dispersion Approximately Balancing Weights: Asymtotic Properties and Practical Considerations," Biometrika, in press.
#' @references Zubizarreta, J. R. (2015), "Stable Weights that Balance Covariates for Estimation with Incomplete Outcome Data," Journal of the American Statistical Association, 110, 910-922.
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
#' bal$bal_std = TRUE
#' 
#' # Solve for the Average Treatment effect among the Treated, ATT (default)
#' bal$bal_alg = FALSE
#' sbwatt.object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal)
#' 
#' # # Solve for a Conditional Average Treatment Effect, CATE
#' # sbwcate.object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "cate", par_tar = "X1 > 1 & X3 <= 0.22"))
#' 
#' # # Solve for the population mean, POP
#' # tar = colMeans(bal_cov)
#' # names(tar) = bal$bal_cov
#' # sbwpop.object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "pop"))
#' 
#' # # Solve for a target population mean, AUX
#' # sbwaux.object = sbw(dat = data_frame, bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "aux", par_tar = tar*1.05))
#' 
#' # # Solve for the ATT using the tuning algorithm
#' # bal$bal_alg = TRUE
#' # bal$bal_sam = 1000
#' # sbwatttun.object = sbw(dat = data_frame, ind = t_ind, out = "Y", bal = bal, 
#' # sol = list(sol_nam = "quadprog"), par = list(par_est = "att", par_tar = NULL))
#' 
#' # Check
#' summarize(sbwatt.object)
#' # summarize(sbwcate.object)
#' # summarize(sbwpop.object)
#' # summarize(sbwaux.object)
#' # summarize(sbwatttun.object)
#' 
#' # Estimate
#' estimate(sbwatt.object)
#' # estimate(sbwcate.object)
#' # estimate(sbwpop.object)
#' # estimate(sbwatttun.object)
#' 
#' # Visualize
#' visualize(sbwatt.object)
#' # visualize(sbwcate.object)
#' # visualize(sbwpop.object)
#' # visualize(sbwaux.object)
#' # visualize(sbwatttun.object)
#' 
#' @export
#'
sbw = function(dat, ind = NULL, out = NULL, bal = list(bal_cov, bal_alg = TRUE, bal_tol, bal_std = TRUE, bal_gri = c(0.0001, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1), bal_sam = 1000), wei = list(wei_sum = TRUE, wei_pos = TRUE), sol = list(sol_nam = "quadprog", sol_dis = FALSE), par = list(par_est = "att", par_tar = NULL)) {
  if (is.null(bal$bal_alg)) {
    bal$bal_alg = TRUE
  } else if (!is.logical(bal$bal_alg)) {
    stop("Please assign a logical to bal$bal_alg.")
  } 
  
  if (is.null(bal$bal_std)) {
    bal$bal_std = TRUE
  } else if (!is.logical(bal$bal_std)) {
    stop("Please assign a logical to bal$bal_std.")
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
  
  if ("weights" %in% colnames(dat)) {
    warning("The 'weights' column in 'dat' will be replaced by the optimal weights. Please change the name of the original 'weights' to avoid unexpectedreplacement.")
  }
  
  if (!par$par_est %in% c("aux")) {
    if (bal$bal_alg == FALSE) {
      object = .sbwfix(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)
    } else if (bal$bal_alg == TRUE) {
      bal$bal_tol = bal$bal_gri
      object = .sbwtun(dat = dat, ind = ind, out = out, bal = bal, wei = wei, sol = sol, par = par)
    }
  } else if (par$par_est %in% c("aux")) {
    if (bal$bal_alg == FALSE) {
      bal$bal_tar = par$par_tar
      object = .sbwauxfix(dat = dat, bal = bal, wei = wei, sol = sol)
    } else if (bal$bal_alg == TRUE) {
      bal$bal_tar = par$par_tar
      object = .sbwauxtun(dat = dat, bal = bal, wei = wei, sol = sol)
    }
    object$ind = ind
    object$out = out
    object$par = par
  }
  return(object)
}