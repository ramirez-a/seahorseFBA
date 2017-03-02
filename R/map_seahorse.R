#'@export map_seahorse
#'@title Map the seahorse fluxes to metabolic reactions
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function takes the matrix of sampled seahorse measurements returned by \code{\link{sample_seahorse}}
#'and maps the fluxes to metabolic reactions.  It returns a matrix with the mapped fluxes of
#'basal oxygen consumption, mitochondrial oxygen consumption, mitochondrial ATP production,
#'mitochondrial ATP leak, and basal extracellular acidification
#'@param model An object of class \code{\link[sybil]{modelorg}}

map_seahorse <- function(x){
  # I may change lactate to protons
  exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "EX_lac_L(e)in", "EX_lac_L(e)ex")  
  
  output_mat <- matrix()
  
  atp_leak <- x[,"ATP_leak"]
  ocr_mito_min <- x[,"OCR_basal"] - x[,"OCR_rotenone"]
  ocr_mito_max <- x[,"OCR_fccp"] - x[,"OCR_rotenone"]
  h_sec <- x[,"PPR_basal"]
}
