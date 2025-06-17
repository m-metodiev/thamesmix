# calculate a matrix used to permute the entries of the parameters
#'
#' Calculate a matrix used to permute the entries of the parameters
#'
#' @param sort_indices      a matrix of indices that will be used for sorting
#' @param num_var_g         number of component mixture parameters
#' @param G                 number of components
#' @importFrom Rfast        rep_row
#' @return a matrix used to permute the entries of the parameters
#'
#' @keywords internal
calc_shift_mat_simple = function(sort_indices,num_var_g,G){
  shift = sort_indices
  for(i in (1:num_var_g)[-num_var_g]){
    shift = cbind(shift,
                  sort_indices +
                    i*Rfast::rep_row(rep(G,ncol(sort_indices)),nrow(shift)))
  }
  return(shift)
}
