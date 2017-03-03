#'@export summarize_seahorse
#'@title Summarize seahorse data by sample and region
#'@author Alfred Ramirez <akram@bu.edu>
#'@description This function takes the long format table exported by the XF Wave software for the 
#'oxidiative stress test and returns a matrix of means and standard deviations for each region
#'of the assay (namely basal, oligomycin, fccp, and rotenone).  The GroupName column must have the sample names.
#'The final measurement of each region is used for the mean and standard deviation across wells.
#'@param x A data.frame
#'@param n The number of regions
# There's definitely a more elegant way of doing this.

summarize_seahorse <- function(x, n=4){
  samples <- unique(x[,"GroupName"])
  
  # I use length(unique(foo)) rather than max(foo) since there may be a time when
  # a whole region should be excluded from the analysis.
  
  m <- length(unique((x[,"Measurement"])))/n #The number measurements in each of the 4 regions
  output_mat <- matrix(nrow=length(samples), ncol=(m*n+1))
  
  for(i in 1:length(samples)){
    mat <- x[x[,"GroupName"] == samples[i],]
    output_mat[i,1] <- samples[i]
    
    
    output_mat[i,m+j] <- mat[1:m, "OCR"]
    output_mat[i,2*m+j] <- mat[1:m, "PPR"]
    output_mat[i,3*m+j] <- mat[1:m, "OCR.Error"]
    output_mat[i,4*m+j] <- mat[1:m, "PPR.Error"]
    
    
    output_mat[i,2] <- mean(mat[1:m,"OCR"]) #mean OCR basal
    output_mat[i,3] <- mean(mat[1:m + m,"OCR"]) #mean OCR oligo
    output_mat[i,4] <- mean(mat[1:m + 2*m,"OCR"]) #mean OCR FCCP
    output_mat[i,5] <- mean(mat[1:m + 3*m,"OCR"]) #mean OCR rotenone
    
    output_mat[i,6] <- mean(mat[1:m,"PPR"]) #mean PPR basal
    output_mat[i,7] <- mean(mat[1:m + m,"PPR"]) #mean PPR oligo
    output_mat[i,8] <- mean(mat[1:m + 2*m,"PPR"]) #mean PPR FCCP
    output_mat[i,9] <- mean(mat[1:m + 3*m,"PPR"]) #mean PPR rotenone
    
    output_mat[i,10] <- sd(mat[1:m,"OCR"]) #sd OCR basal
    output_mat[i,11] <- sd(mat[1:m + m,"OCR"]) #sd OCR oligo
    output_mat[i,12] <- sd(mat[1:m + 2*m,"OCR"]) #sd OCR FCCP
    output_mat[i,13] <- sd(mat[1:m + 3*m,"OCR"]) #sd OCR rotenone
    
    output_mat[i,14] <- sd(mat[1:m,"PPR"]) #sd PPR basal
    output_mat[i,15] <- sd(mat[1:m + m,"PPR"]) #sd PPR oligo
    output_mat[i,16] <- sd(mat[1:m + 2*m,"PPR"]) #sd PPR FCCP
    output_mat[i,17] <- sd(mat[1:m + 3*m,"PPR"]) #sd PPR rotenone
  }
  colnames(output_mat) <- c("Sample", "OCR_basal","OCR_oligo", "OCR_fccp", "OCR_rotenone",
                            "PPR_basal", "PPR_oligo", "PPR_fccp", "PPR_rotenone",
                            "OCR_basal_sd","OCR_oligo_sd", "OCR_fccp_sd", "OCR_rotenone_sd", 
                            "PPR_basal_sd", "PPR_oligo_sd", "PPR_fccp_sd", "PPR_rotenone_sd")
  output_mat
}
