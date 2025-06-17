#' Extends the parameter to include the last proportion
#'
#' The last proportion parameters is redundant, since it is equal to 1 minus
#' the sum of all other proportion parameters. This function ats the last
#' proportion parameter back to the parameter matrix.
#'
#' @param param_test  a matrix of parameters
#' @param G           number of components
#' @return the matrix param_test with one additional column
#'
#' @keywords internal
extend_param = function(param_test, G){
  if((ncol(param_test)/G) > 1){
    if(G==2){
      param_test_extended = cbind(param_test,1-param_test[,ncol(param_test)])
    } else{
      param_test_extended =
        cbind(param_test,
              1-rowSums(param_test[,(ncol(param_test)-(G-2)):ncol(param_test)]))
    }
  } else{
    param_test_extended = param_test
  }
  return(param_test_extended)
}
