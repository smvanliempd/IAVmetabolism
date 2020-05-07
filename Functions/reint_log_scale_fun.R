reint.log.scale.fun <- function( reint.adjust.out ) {
  
  dat <- reint.adjust.out$data
  
  # log trans all
  dat[ !is.na(Area_reint_mfc_qc) , reint_log := log10(Area_reint_mfc_qc) ]
  
  # standardize 
  dat[ Sample.Group != "QC" , reint_log_mean := mean(reint_log, na.rm = T), by = Feature ]
  dat[ Sample.Group != "QC" , reint_log_sd   := sd(reint_log, na.rm = T)    , by = Feature]
  dat[ Sample.Group != "QC" , reint_log_stand:= (reint_log - reint_log_mean)/reint_log_sd ]
  
  # out
  reint.adjust.out$data <- dat
  
  return(reint.adjust.out)
  
}