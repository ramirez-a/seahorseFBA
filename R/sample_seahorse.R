#'@title Sample from the seahorse data
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function takes the data.frame returned by the function \code{\link{summarize_seahorse}} and returns
#'A data.frame with sampled measurements
#'@param x A data.frame

sample_seahorse <- function(x){
  bsamples <- length(unique(x$Sample))
  output_mat <- cbind(x, matrix(0, ncol=2, nrow=bsamples))
  colnames(output_mat)[(ncol(output_mat)-1):ncol(output_mat)] <- c("ATP_flux", "ATP_leak")
  
  for(i in 1:bsamples){
    output_mat[i,"OCR_basal"] <- rnorm(1,mean=x[i,"OCR_basal"],sd=x[i,"OCR_basal_sd"])
    output_mat[i,"OCR_oligo"] <- rnorm(1,mean=x[i,"OCR_oligo"],sd=x[i,"OCR_oligo_sd"])
    output_mat[i,"OCR_fccp"] <- rnorm(1,mean=x[i,"OCR_fccp"],sd=x[i,"OCR_fccp_sd"])
    output_mat[i,"OCR_rotenone"] <- rnorm(1,mean=x[i,"OCR_rotenone"],sd=x[i,"OCR_rotenone_sd"])
    
    output_mat[i,"PPR_basal"] <- rnorm(1,mean=x[i,"PPR_basal"],sd=x[i,"PPR_basal_sd"])
    output_mat[i,"PPR_oligo"] <- rnorm(1,mean=x[i,"PPR_oligo"],sd=x[i,"PPR_oligo_sd"])
    output_mat[i,"PPR_fccp"] <- rnorm(1,mean=x[i,"PPR_fccp"],sd=x[i,"PPR_fccp_sd"])
    output_mat[i,"PPR_rotenone"] <- rnorm(1,mean=o[i,"PPR_rotenone"],sd=x[i,"PPR_rotenone_sd"])
    
    output_mat[i,"ATP_flux"] <- 4.6*(output_mat[i,"OCR_basal"] - output_mat[i,"OCR_oligo"])
    output_mat[i,"ATP_leak"] <- 4.6*(output_mat[i,"OCR_oligo"] - output_mat[i,"OCR_rotenone"])
  }  
  output_mat
}
