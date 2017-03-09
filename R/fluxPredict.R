#'@export fluxPredict
#'@title Make predictions from seahorse data by sampling and minimizing total flux
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function integrates the sampled seahorse measurements as constraints into the specified model,
#'minimizes total flux for each sample, and returns a matrix of reactions x samples where the entries are predicted fluxes.
#'@param model An object of class \code{\link[sybil]{modelorg}}
#'@param seahorse_data A data.frame returned by \code{\link{map_seahorse}}
#'@param biomass_est Estimated biomass flux.  Not required, by may be useful.

fluxPredict <-function(model, seahorse_data, biomass_est=0){
  #Sample
  library(foreach)
  library(doMC)
  registerDoMC(cores=cores)
  
  ba_output <- foreach(i=1:nsamples, .combine=cbind) %dopar% {
    tissue <- "hBA"
    atps4m_ba <- -1
    atp_leak_ba <- -1
    lpsolution <- 0
    
    #Before MTF, we need to make sure the sampled inputs and FBA make sense
    while(atps4m_ba <= 0 | atp_leak_ba <= 0 | lpsolution != 1) {
      adipocyte_data <- sample_seahorse_tseng(adipocyte_data_original)
      ocr_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"]
      
      #atps4m_ba must be postive
      atps4m_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_flux"] 
      
      atp_leak_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"ATP_leak"]
      ocr_mito_min_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_basal"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
      ocr_mito_max_ba <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_fccp"] - adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"OCR_rotenone"]
      h_sec <- adipocyte_data[adipocyte_data[,"Tissue"] == tissue,"PPR_basal"]
      biomass_ba <- log(2)/76.8
      
      ba_lb <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_min_ba, biomass_ba,0, h_sec)
      ba_ub <- c(ocr_ba,0, atps4m_ba, atp_leak_ba, ocr_mito_max_ba, biomass_ba,0, h_sec)
      ba_model <- changeBounds(pread_model, react= exp_coefs, lb=ba_lb, ub=ba_ub)
      ba_model <- changeObjFunc(ba_model, "biomass_reaction", obj_coef=1)
      
      ba_test <- changeObjFunc(ba_model, "biomass_reaction", obj_coef=1)
      ba_test <- changeBounds(ba_test, c("PALFATPtc", "RTOTAL2FATPc_pmt", "RTOTALFATPc_pmt", "RTOTAL3FATPc_pmt"), lb=rep(0,4), ub=rep(0,4))
      ba_test_flux <- optimizeProb(ba_test, algorithm="fba", lpdir="max")
      
      lpsolution <- ba_test_flux@lp_stat
    }
    
    ba_fluxes <- optimizeProb(ba_model, alg="mtf", mtfobj=mod_obj(ba_test_flux))
    ba_fluxes_df <- getFluxDist(ba_fluxes)
    ba_output <- ba_fluxes_df
  }
  rownames(ba_output) <- secretory_reactions
  
}