# Median fold change adjustment
mfc.fun    <- function( clean.output ) {
  
  # takes data from the clean.fun() output
  
  ################################# median fold change (MFC) normalization #####################################
  # In order to make proper comparisons between samples, feature signals should be normalized between samples. #
  # Therefore we use the *median fold-change* method (MFC normalization). In short, a fold-change (FC) for     #
  # features per sample is calculated with respect to a reference sample. The reference sample is chosen to      #
  # be the one with the least missing values. Then, per sample, the median is taken over all the FC features.  #
  # This value is used to normalize all the features within a sample by deviding it by the Feature signal.     #    
  ##############################################################################################################
  
  # feature signal data
  dat <- clean.output$data
  
  ## MFC calculations
  # Select sample with least NA's, QCs excluded
  mfc.ref <- dat[Sample.Group != "QC", sum(!is.na(Signal) ), by =  File.Name ][V1 == max(V1), File.Name] 
  
  # Choose 1st sample in case of ties
  mfc.ref <- mfc.ref[1]                                                                                  
  
  # Calulate fold change per feature
  dat[ , fold_change := Signal/Signal[ File.Name == mfc.ref], by = Feature ]                              
  
  # Calculte median fold change  for each sample
  dat[ ,  mfc_factor := median(fold_change, na.rm = T), by = File.Name]                                  
  
  # Devide signal by mfc factor
  dat[ ,  Signal_mfc := Signal/mfc_factor]                                                              
  
  # out
  clean.output$data <- dat
  clean.output$mfc_ref <- mfc.ref
  
  return(clean.output)
  
}