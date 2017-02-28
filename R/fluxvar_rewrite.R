#'@title Rewrite fluxVar result into a dataframe
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function rewrites an object of \code{\link[sybil]{optsol}}
#'@param model A object of class \code{\link[sybil]{modelorg}}
#'@param reactions A character vector of reaction ids
#'@param solver The solver to use.  Default SYBIL_SETTINGS("SOLVER")

#This function takes the results from fluxVar and rewrites them into a convenient matrix
fluxvar_rewrite <- function(o){
  results_vector <- o@lp_obj
  reactions_tested <- o@react@react_id
  output_mat <- matrix(results_vector, byrow=F, ncol=2, nrow=length(reactions_tested))
  colnames(output_mat) <- c("Min", "Max")
  rownames(output_mat) <- reactions_tested
  output_mat
}
