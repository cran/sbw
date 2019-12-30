# functions for solving linear programming
.as.sparseMatrix <- function(slam_matrix) {
  retval <-  Matrix::sparseMatrix(i=as.numeric(slam_matrix$i),
                                  j=as.numeric(slam_matrix$j),
                                  x=as.numeric(as.character(slam_matrix$v)),
                                  dims=c(slam_matrix$nrow, 
                                         slam_matrix$ncol),
                                  dimnames = dimnames(slam_matrix),
                                  giveCsparse = TRUE)
}

# solver: pogs
.sbwpripogs = function(problemparameters.object, params) {
  normalize = problemparameters.object$normalize
  nor = problemparameters.object$nor
  n = problemparameters.object$n
  Amat = problemparameters.object$Amat
  bvec = problemparameters.object$bvec
  cvec = problemparameters.object$cvec
  lb = problemparameters.object$lb
  ub = problemparameters.object$ub
  sense = problemparameters.object$sense
  rm(problemparameters.object)
  Amat = as.matrix(Amat)
  
  if (nor == "l_1" | nor == "l_inf") {
    stop("Minimizing the l_1 or l_inf norm is not a valid option with pogs. Please use one of cplex, gurobi and mosek.")
  }
  if (nor == "l_2") {
    # The primal problem is
    # minimize: c'x + 0.5x'Qx
    # subject to:
    # 	A x <= b
    #   (and sum(x) = 1 if normalize = 1)
    # with bounds:
    # 	lb <= x <= ub
    # The solver pogs has the form
    # minimize: f(y) + g(x), s.t. y = Ax,
    # where f and g are convex, separable, and take the form
    #   c h(a x - b) + d x + e x^2
    #   where a, b and d are real, c and d are non-negative and h is one of 16 convex functions
    if(normalize == 1) {
      meq = sum(sense == "E")
      nbL = nrow(Amat) - meq
      nbC = ncol(Amat)
      # for balance conditions
      f = list(h = c(pogs::kIndLe0(nbL), pogs::kIndEq0(meq)), b = bvec, a = 1, c = 1, d = 0, e = 0)
      # for positive weights and l2 norm penalty
      g = list(h = pogs::kIndGe0(), c = 1, a = 1, d = 0, e = 1)
      cat(format("  POGS optimizer is opening..."), "\n")
      cat(format("  Finding the optimal weights..."), "\n")
      out = pogs::pogs(Amat, f, g, params = params)
    } else {
      nbL = nrow(Amat)
      nbC = ncol(Amat)
      Amat2 = matrix(0, (nbL+1), (nbC+1))
      Amat2[1:nbL, 1:nbC] = Amat
      Amat2[(nbL+1),] = c(rep(1,nbC), -1)
      f = list(h = c(pogs::kIndLe0(nbL), pogs::kIndEq0(1)), b = c(bvec,0), a = 1, c = 1, d = 0, e = 0)
      if (lb[1] == 0) {
        g = list(h = pogs::kIndGe0(), c = 1, a = 1, d = 0, e = c(rep(1, nbC), -1/n))
      }
      if (lb[1] == -Inf) {
        g = list(h = pogs::kZero(), c = 1, a = 1, d = 0, e = c(rep(1, nbC), -1/n))
      }
      
      # nbL = nrow(Amat)
      # nbC = ncol(Amat)
      # # for balance conditions
      # f = list(h = c(pogs::kIndLe0(nbL)), b = bvec, a = 1, c = 1, d = 0, e = 0)
      # # for positive weights and variance
      # if (lb[1] == 0) {
      #   g = list(h = pogs::kIndGe0(), c = 1, a = 1, d = 0, e = 1)
      # }
      # # for positive weights and variance
      # if (lb[1] == -Inf) {
      #   g = list(h = pogs::kZero(), c = 1, a = 1, d = 0, e = 1)
      # }
      cat(format("  POGS optimizer is opening..."), "\n")
      cat(format("  Finding the optimal weights..."), "\n")
      out = pogs::pogs(Amat2, f, g, params = params)
    }
  }
  # Get status code
  status = out$status
  if (status == 0) {
    status = "optimal"
    cat(format("  Optimal weights found."), "\n")
    # Get weights
    weights = (out$x)[1:n]
    if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
      weights[weights < 0] = 0
    }
    if (normalize == 1) {
      weights = weights/sum(weights)
    }
    # Get objective
    obj_total = 2*out$optval - n*sum(weights/n)^2
    # Get dual table
    dual_table = NULL
  } else {
    cat(format("  Optimal weights not found."), "\n")
    stop(paste("pogs status code = ", status))
  }
  
  return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
}

# solver: cplex
.sbwpricplex = function(problemparameters.object, trace) {
  normalize = problemparameters.object$normalize
  nor = problemparameters.object$nor
  n = problemparameters.object$n
  Amat = problemparameters.object$Amat
  bvec = problemparameters.object$bvec
  cvec = problemparameters.object$cvec
  lb = problemparameters.object$lb
  ub = problemparameters.object$ub
  sense = problemparameters.object$sense
  var_type = problemparameters.object$var_type
  rm(problemparameters.object)
  
  cat(format("  CPLEX optimizer is opening..."), "\n")
  cat(format("  Finding the optimal weights..."), "\n")
  if (nor == "l_1" | nor == "l_inf") {
    out = Rcplex::Rcplex(cvec, Amat, bvec, lb = lb, ub = ub, sense = sense, vtype = var_type, n = 1, control = list(trace = trace))
  }
  if (nor == "l_2") {
    # The primal problem is
    # minimize: c'x + 0.5x'Qx
    # subject to:
    # 	A x <= b
    #   (and sum(x) = 1 if normalize = 1)
    # with bounds:
    # 	lb <= x <= ub
    # In the solver cplex, 
    #   cvec := c
    #   Amat := A
    #   bvec := b
    #   Qmat := Q
    #   lb <= x <= ub
    # For l_2 norm,
    #   c = 0
    #   Q = 2*(diag(n) - 1/n)
    Qmat = diag(n) - 1/n
    out = Rcplex::Rcplex(cvec = cvec, Amat = Amat, bvec = bvec, Qmat = 2*Qmat, lb = lb, ub = ub, sense = sense, vtype = var_type, n = 1, control = list(trace = trace))
  }
  # Get status code
  status = out$status
  if (status == 1) {
    status = "optimal"
    cat(format("  Optimal weights found."), "\n")
  } else if (status %in% c(2, 3, 4)) {
    # status = "infeasible or unbounded"
    cat(format("  Problem ill-posed."), "\n")
    message("  Please find the cplex status code on https://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.3/ilog.odms.cplex.help/refcallablelibrary/macros/homepagesolutionstatus.html")
    stop(paste("CPLEX status code = ", status))
  } else {
    cat(format("  Optimal weights not found."), "\n")
    message("  Please find the cplex status code on https://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.3/ilog.odms.cplex.help/refcallablelibrary/macros/homepagesolutionstatus.html")
    stop(paste("CPLEX status code = ", status))
  }
  # Get weights
  weights = (out$xopt)[1:n]
  if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
    weights[weights < 0] = 0
  }
  if (normalize == 1) {
    weights = weights/sum(weights)
  }
  # Get objective
  obj_total = out$obj
  # Get dual table
  dual_vars = out$extra$lambda[names(bvec) != ""]
  dual_names = names(bvec)[names(bvec) != ""]
  if (!is.null(dual_vars)) {
    names(dual_vars) = dual_names
    dual_table  = .dualtable(dual_vars)
  } else dual_table = NULL
  
  return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
}  
   
# solver: quadprog
.sbwpriquadprog = function(problemparameters.object) {
  normalize = problemparameters.object$normalize
  nor = problemparameters.object$nor
  n = problemparameters.object$n
  Amat = problemparameters.object$Amat
  bvec = problemparameters.object$bvec
  cvec = problemparameters.object$cvec
  lb = problemparameters.object$lb
  ub = problemparameters.object$ub
  sense = problemparameters.object$sense
  rm(problemparameters.object)
  
  if (nor == "l_1" | nor == "l_inf") {
    stop("Minimizing the l_1 or l_inf norm is not a valid option with quadprog. Please use one of cplex, gurobi and mosek.")
  }
  if (nor == "l_2") {
    # The primal problem is
    # minimize: c'x + 0.5x'Qx
    # subject to:
    # 	A x <= b
    #   (and sum(x) = 1 if normalize = 1)
    # with bounds:
    # 	lb <= x <= ub
    # In the solver quadprog, 
    #   dvec := -c
    #   Amat := -(A, -I)'
    #   bvec := (-b, 0)
    #   Dmat := Q
    # The first meq constraints are treated as equality constraints
    #   if normalize = 1, meq = 1
    #   Amat = (1, Amat)
    #   bvec = (1, bvec)
    # For l_2 norm,
    #   c = 0
    #   Q = 2*diag(n)
    Dmat = 2*diag(n)
    Amat = as.matrix(Amat)

    if (normalize == 1) {
      # Intervert the first line with the last line of Amat and bvec
      meq = sum(sense == "E")
      tmp = bvec[(length(bvec) - meq + 1):length(bvec)]
      bvec[(meq + 1):length(bvec)] = bvec[1:(length(bvec) - meq)]
      bvec[1:meq] = tmp
      tmp = Amat[(nrow(Amat) - meq + 1):nrow(Amat),]
      Amat[(meq + 1):nrow(Amat),] = Amat[1:(nrow(Amat) - meq),]
      Amat[1:meq,] = tmp
    }
    
    if (lb[1] == -Inf) {
      Amat2 = Amat
      bvec2 = bvec
    }
    if (lb[1] == 0) {
      # Add the non-negativity constraints
      Amat2 = matrix(0, nrow(Amat) + ncol(Amat), ncol(Amat))
      Amat2[1:nrow(Amat), 1:ncol(Amat)] = Amat
      Amat2[(1 + nrow(Amat)):(nrow(Amat) + ncol(Amat)), 1:ncol(Amat)] = -1*diag(ncol(Amat))
      bvec2 = rep(0, (length(bvec) + n))
      bvec2[1:length(bvec)] = bvec
    }
    # Adapt matrices to quadprog
    # minimize: -dvec'x + 0.5x'Dmat x
    # subject to:
    #   Amat x >= bvec
    # the first meq constraints are treated as equality constraints, meq
    Amat4 = -t(Amat2)
    bvec4 = -bvec2

    cat(format("  quadprog optimizer is opening..."), "\n")
    cat(format("  Finding the optimal weights..."), "\n")
    out = quadprog::solve.QP(Dmat, dvec = cvec, Amat4, bvec4, meq = 1)
  }
  # Get weights
  weights = (out$solution)[1:n]
  if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
    weights[weights < 0] = 0
  }
  if (normalize == 1) {
    weights = weights/sum(weights)
  }
  # Get objective
  obj_total = out$value - n*mean(weights)^2
  if (is.na(obj_total) | is.nan(obj_total)| is.null(obj_total)) {
    cat(format("  Optimal weights not found."), "\n")
  } else cat(format("  Optimal weights found."), "\n")
  # No status code
  status = NA
  # Get dual table
  if (normalize == 1) {
    dual_vars = (-out$Lagrangian[(meq+1):(nrow(Amat))])
  } else dual_vars = (-out$Lagrangian[1:(nrow(Amat))])
  dual_names = names(bvec)[names(bvec) != ""]
  if (!is.null(dual_vars)) {
    names(dual_vars) = dual_names
    dual_table  = .dualtable(dual_vars)
  } else dual_table = NULL
  
  return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
}  


# # solver: quadprogpp
# .sbwpriquadprogpp = function(problemparameters.object) {
  # normalize = problemparameters.object$normalize
  # nor = problemparameters.object$nor
  # n = problemparameters.object$n
  # Amat = problemparameters.object$Amat
  # bvec = problemparameters.object$bvec
  # cvec = problemparameters.object$cvec
  # lb = problemparameters.object$lb
  # ub = problemparameters.object$ub
  # sense = problemparameters.object$sense
  # rm(problemparameters.object)
  
  # if (nor == "l_1" | nor == "l_inf") {
    # stop("Minimizing the l_1 or l_inf norm is not a valid option with quadprogpp. Please use one of cplex, gurobi and mosek.")
  # }
  # if (nor == "l_2") {
    # # The primal problem is
    # # minimize: c'x + 0.5x'Qx
    # # subject to:
    # # 	A x <= b
    # #   (and sum(x) = 1 if normalize = 1)
    # # with bounds:
    # # 	lb <= x <= ub
    # # In the solver quadprogpp, 
    # #   G := Q
    # #   g0 := c
    # #   CI := -(A', -I)
    # #   ci0 := (b, 0)
    # #   if normalize = 1, meq = 1
    # #   CE = matrix(1, n, 1)
    # #   ce0 = c(-1)

    # G = 2*diag(n)
    # g0 = cvec
    # Amat = as.matrix(Amat)
    
    # if (normalize == 1) {
      # meq = sum(sense == "E")
      # ci0 = bvec[1:(length(bvec) - meq)]
      # CI = -t(Amat[1:(nrow(Amat) - meq),])
      # ce0 = c(-1)
      # CE = matrix(1, n, 1)
    # }
    
    # if (lb[1] == 0) {
      # # Add the non-negativity constraints
      # CI = cbind(CI, diag(n))
      # ci0 = c(ci0, rep(0, n))
    # }
    
    # cat(format("  quadprogpp optimizer is opening..."), "\n")
    # cat(format("  Finding the optimal weights..."), "\n")
    # out = quadprogpp::QP.Solve(G, g0, CI, ci0, CE, ce0)
  # }
  # # Get weights
  # weights = out[1:n]
  # if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
    # weights[weights < 0] = 0
  # }
  # if (normalize == 1) {
    # weights = weights/sum(weights)
  # }
  # # Get objective
  # obj_total = NA
  # if (is.na(weights[1])) {
    # cat(format("  Optimal weights not found."), "\n")
  # } else cat(format("  Optimal weights found."), "\n")
  # # No status code
  # status = NA
  # # Get dual table
  # dual_table = NULL
  
  # return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
# }  


# solver: gurobi
.sbwprigurobi = function(problemparameters.object, params) {
  normalize = problemparameters.object$normalize
  nor = problemparameters.object$nor
  n = problemparameters.object$n
  Amat = problemparameters.object$Amat
  bvec = problemparameters.object$bvec
  cvec = problemparameters.object$cvec
  lb = problemparameters.object$lb
  ub = problemparameters.object$ub
  sense = problemparameters.object$sense
  var_type = problemparameters.object$var_type
  rm(problemparameters.object)
  
  # The primal problem is
  # minimize: c'x + 0.5x'Qx
  # subject to:
  # 	A x <= b
  #   (and sum(x) = 1 if normalize = 1)
  # with bounds:
  # 	lb <= x <= ub
  # In solver gurobi, 
  #   obj := c
  #   A := A
  #   rhs := b
  #   Q := Q
  model = list()
  model$modelsense = 'min'
  model$obj = cvec
  model$A = Amat
  model$sense = rep(NA, length(sense))
  model$sense[sense=="E"] = '='
  model$sense[sense=="L"] = '<='
  model$sense[sense=="G"] = '>='
  model$rhs = bvec
  model$vtype = var_type
  model$lb = lb
  model$ub = ub
    
  cat(format("  Gurobi optimizer is opening..."), "\n")
  cat(format("  Finding the optimal weights..."), "\n")
  if (nor == "l_1" | nor == "l_inf") {
    out = gurobi::gurobi(model, params)
  }
  if (nor == "l_2") {
    model$Q = diag(n) - 1/n
    out = gurobi::gurobi(model, params)
  }
  # Get status code
  status = out$status
  if (status == "OPTIMAL") {
    status = "optimal"
    cat(format("  Optimal weights found."), "\n")
  } else {
    cat(format("  Optimal weights not found."), "\n")
    message("  Please find the gurobi status code on http://www.gurobi.com/documentation/8.0/refman/optimization_status_codes.html")
    stop(paste("Gurobi status code = ", status))
  }
  # Get weights
  weights = (out$x)[1:n]
  if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
    weights[weights < 0] = 0
  }
  if (normalize == 1) {
    weights = weights/sum(weights)
  }
  # Get objective
  obj_total = out$obj
  # Get dual table
  dual_vars = out$pi[names(bvec) != ""]
  dual_names = names(bvec)[names(bvec) != ""]
  if (!is.null(dual_vars)) {
    names(dual_vars) = dual_names
    dual_table  = .dualtable(dual_vars)
  } else dual_table = NULL
  
  return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
}

# solver = mosek
.sbwprimosek = function(problemparameters.object, verbose) {
  normalize = problemparameters.object$normalize
  nor = problemparameters.object$nor
  n = problemparameters.object$n
  Amat = problemparameters.object$Amat
  bvec = problemparameters.object$bvec
  cvec = problemparameters.object$cvec
  lb = problemparameters.object$lb
  ub = problemparameters.object$ub
  sense = problemparameters.object$sense
  
  rm(problemparameters.object)
  
  if (nor == "l_2") {
    # The primal problem is
    # minimize: c'x + 0.5x'Qx
    # subject to:
    # 	A x <= b
    #   (and sum(x) = 1 if normalize = 1)
    # with bounds:
    # 	lb <= x <= ub
    # Solver mosek has the form, 
    # minimize: f'x + 0.5x'(F'F)x
    # subject to:
    # 	A x <= b
    #  	Aeq x = beq
    # with bounds:
    # 	lb <= x <= ub
    # For l_2 norm,
    #   f = 0
    #   F = sqrt(2)*diag(n)
    mos = list()
    if (normalize == 1) {
      meq = sum(sense == "E")
      mos$Aeq = as.matrix(Amat[(nrow(Amat) - meq + 1):nrow(Amat),], nrow = meq)
      mos$beq = bvec[(length(bvec) - meq + 1):length(bvec)]
      # Amat = as.matrix(Amat)
      Amat = .as.sparseMatrix(Amat)
      mos$A = Amat[-((nrow(Amat) - meq + 1):nrow(Amat)),]
      mos$b = bvec[-((length(bvec) - meq + 1):length(bvec))]
    } else {
        # Amat = as.matrix(Amat)
        Amat = .as.sparseMatrix(Amat)
        mos$Aeq = NA
        mos$beq = NA
        mos$A = Amat
        mos$b = bvec
    }
    mos$F = Matrix::Diagonal(n)*sqrt(2)
    mos$f = cvec
    mos$lb = lb
    mos$ub = ub
      
    prob = Rmosek::mosek_qptoprob(F = mos$F, f = mos$f, A = mos$A, b = mos$b, Aeq = mos$Aeq, beq = mos$beq, lb = mos$lb, ub = mos$ub)
    cat(format("  Mosek optimizer is opening..."), "\n")
    cat(format("  Finding the optimal weights..."), "\n")
    out = Rmosek::mosek(prob, opts = list(verbose = verbose))
  }
  if (nor == "l_1" | nor == "l_inf") {
    # For l_1 and l_inf norm
    # minimize: f'x
    #   subject to:
    #   A x <= b
    #   Aeq x = beq
    # with bounds:
    #  lb <= x <= ub
      
    mos = list()
    if (normalize == 1) {
      meq = sum(sense == "E")
      mos$Aeq = as.matrix(Amat[(nrow(Amat)-meq+1):nrow(Amat),], nrow = meq)
      mos$beq = bvec[(length(bvec)-meq+1):length(bvec)]
      Amat = .as.sparseMatrix(Amat)
      mos$A = Amat[-((nrow(Amat)-meq+1):nrow(Amat)),]
      mos$b = bvec[-((length(bvec)-meq+1):length(bvec))]
    } else {
      mos$Aeq = NA
      mos$beq = NA
      mos$A = Amat
      mos$b = bvec
    }
    mos$f = cvec
    mos$lb = lb
    mos$ub = ub
      
    prob = Rmosek::mosek_lptoprob(f = mos$f, A = mos$A, b = mos$b, Aeq = mos$Aeq, beq = mos$beq, lb = mos$lb, ub = mos$ub)
    cat(format("  Mosek optimizer is opening..."), "\n")
    cat(format("  Finding the optimal weights..."), "\n")
    out = Rmosek::mosek(prob, opts = list(verbose = verbose))
  }
  # Get status
  code = out$response$code
  solsta = out$sol$it$solsta
  if (!is.null(solsta)) {
    if (code == 0) {
      # Get weights
      weights = out$sol$itr$xx[1:n]
      if (solsta == "OPTIMAL") {
        status = "optimal"
        cat(format("  Optimal weights found."), "\n")
        obj_total = as.numeric(prob$c %*% out$sol$itr$xx - n*mean(weights)^2)
        if (lb[1] == 0 & min(weights, na.rm = TRUE) < 0) {
          weights[weights < 0] = 0
        }
        if (normalize == 1) {
          weights = weights/sum(weights)
        }
        # Dual variable for upper constraints
        dual_vars = (-out$sol$itr$suc[1:length(bvec)])[names(bvec) != ""]
        dual_names = names(bvec)[names(bvec) != ""]
        if (!is.null(dual_vars)) {
          names(dual_vars) = dual_names
          dual_table  = .dualtable(dual_vars)
        } else dual_table = NULL
      } else {
        status = solsta
        cat(format("  Problem ill-posed."), "\n")
        message("  Please find the mosek status code on https://docs.mosek.com/8.1/capi/response-codes.html")
        stop(paste("Mosek status = ", status))
      }
    } else {
      cat(format("  Optimal weights not found."), "\n")
      message("  Please find the mosek status code on https://docs.mosek.com/8.1/capi/response-codes.html")
      stop(paste("Mosek status code = ", code))
      } 
  } else {
    cat(format("  Optimal weights not found."), "\n")
    message("  Please find the mosek status code on https://docs.mosek.com/8.1/capi/response-codes.html")
    stop(paste("Mosek status code = ", code))
  }  
    
  return(list(obj_total = obj_total, status = status, weights = weights, dual_table = dual_table))
}