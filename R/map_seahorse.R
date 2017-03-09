#'@export map_seahorse
#'@title Map the seahorse fluxes to metabolic reactions
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function takes the matrix of sampled seahorse measurements returned by \code{\link{sample_seahorse}}
#'and maps the fluxes to metabolic reactions. It returns a matrix with the mapped fluxes of
#'basal oxygen consumption, mitochondrial oxygen consumption, mitochondrial ATP production,
#'mitochondrial ATP leak, and basal extracellular acidification
#'@param model An object of class \code{\link[sybil]{modelorg}}

map_seahorse <- function(x){
  # These reaction names are specific to Recon 2.1x.  Probably won't work for Recon 2.2.
  exp_coefs <- c("EX_o2(e)in","EX_o2(e)ex","ATPS4m", "DM_atp_m_", "O2tm", "O2tm_ub", "EX_lac_L(e)in", "EX_lac_L(e)ex")  
  
  output_mat <- matrix(0, ncol=ncol(x), nrow=length(exp_coefs))
  colnames(output_mat) <- colnames(x)
  rownames(output_mat) <- exp_coefs
  
  output_mat["EX_o2(e)in",] <- x["OCR_basal",]
  output_mat["EX_o2(e)ex",] <- 0
  output_mat["ATPS4m",] <- 4.6*(x["OCR_basal",] - x["OCR_oligo",])
  output_mat["DM_atp_m_",] <- 4.6*(x["OCR_oligo",] - x["OCR_rotenone",])
  output_mat["O2tm",] <- x["OCR_basal",] - x["OCR_rotenone",]
  output_mat["O2tm_ub",] <- x["OCR_fccp",] - x["OCR_rotenone",]
  
  #Lactate and extracellular acid may be negative and the associated reactions have different 
  #stoichiometries due to the total carbon constraint, thus a if statement is needed to check to determine which one to set.
  #Ideally, there should never be a case where the PPR negative (except a few edges cases in biology).
  for(i in 1:ncol(x)){
    if(x["PPR_basal",i] < 0){
      output_mat["EX_lac_L(e)in",] <- abs(x["PPR_basal",i])
      output_mat["EX_lac_L(e)ex",] <- 0
    } else{
      output_mat["EX_lac_L(e)in",] <- 0
      output_mat["EX_lac_L(e)ex",] <- x["PPR_basal",]      
    }
  }

  if(any(output_mat < 0 )){
    warning("Negative values were found. Some samples may have values inconsistent with biological expectation")
  }
  
  output_mat
}
