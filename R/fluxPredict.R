#'@export fluxPredict
#'@title Make predictions from seahorse data by sampling and minimizing total flux
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function integrates the seahorse measurements as constraints into the specified model,
#'minimizes total flux for each sample, and returns a matrix of reactions x samples where the entries are predicted fluxes.
#'@param model An object of class \code{\link[sybil]{modelorg}}
#'@param seahorse_data A data.frame
#'@param nsamples The number of samples to perform

fluxPredict <-function(model, seahorse_data, nsamples=150){
  
}