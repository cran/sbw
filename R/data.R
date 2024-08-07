#' The Lalonde data set
#' 
#' Data set from the National Supported Work Demonstration 
#' (Lalonde 1986, Dehejia and Wahba 1999). 
#' This data set is publicly available at 
#' \url{https://users.nber.org/~rdehejia/data/.nswdata2.html}.
#' 
#' @usage data(lalonde)
#' 
#' @format A data frame with 614 observations, corresponding to 185 treated 
#' and 429 control subjects, and 10 variables.  
#' The treatment assignment indicator is the first variable of the data frame; 
#' the next eight columns are the covariates; the last column is the outcome:
#' \describe{
#'       \item{treatment}{the treatment assignment indicator (1 if treated, 0 otherwise)}
#'       \item{age}{a covariate, measured in years}
#'       \item{education}{a covariate, measured in years}
#'       \item{black}{a covariate indicating race (1 if black, 0 otherwise)}
#'       \item{hispanic}{a covariate indicating race (1 if Hispanic, 0 otherwise)}
#'       \item{married}{a covariate indicating marital status (1 if married, 0 otherwise)}
#'       \item{nodegree}{a covariate indicating high school diploma (1 if no degree, 0 otherwise)}
#'       \item{re74}{a covariate, real earnings in 1974}
#'       \item{re75}{a covariate, real earnings in 1975}
#'       \item{re78}{the outcome, real earnings in 1978}
#'       }
#' @source{\url{https://users.nber.org/~rdehejia/data/.nswdata2.html}} 
#' @references Dehejia, R., and Wahba, S. (1999), "Causal Effects in Nonexperimental Studies: Reevaluating the Evaluation of Training Programs," \emph{Journal of the American Statistical Association}, 94, 1053-1062. 
#' @references Lalonde, R. (1986), "Evaluating the Econometric Evaluations of Training Programs," \emph{American Economic Review}, 76, 604-620.
"lalonde"
