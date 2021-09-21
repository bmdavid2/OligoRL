
add_objective_constraint <- function(model, result, frac=1.0, clear_obj=TRUE) {
  m <- nrow(model$A)
  A <- matrix(0, nrow=m+1, ncol=ncol(model$A))
  A[1:m,1:ncol(model$A)] <- model$A
  A[m+1, ] <- model$obj
  model$A <- A
  model$sense[m+1] <- ">"
  model$rhs[m+1] <- frac * result$objval

  if (clear_obj) {
    model$obj[] <- 0
  }

  return(model)
}

#' Test the Gurobi MILP solver
#'
#' @export
test_gurobi <- function() {
  model <- list()

  model$A          <- matrix(c(1,2,3,1,1,0), nrow=2, ncol=3, byrow=T)
  model$obj        <- c(1,1,2)
  model$modelsense <- "max"
  model$rhs        <- c(4,1)
  model$sense      <- c('<', '>')
  model$vtype      <- 'B'

  params <- list(OutputFlag=1)

  result <- gurobi::gurobi(model, params)
}

