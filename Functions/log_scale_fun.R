# log transform and standardize signals for modeling
log.scale.fun <- function( fill.out ) {
  
  dat <- fill.out$data
  
  # log trans all
  dat[ !is.na(Signal_mfc_qc) , log_signal := log10(Signal_mfc_qc) ]
  
  # standardize 
  dat[ Sample.Group != "QC" , log_signal_mean := mean(log_signal, na.rm = T), by = Feature ]
  dat[ Sample.Group != "QC" , log_signal_sd   := sd(log_signal, na.rm = T)    , by = Feature]
  dat[ Sample.Group != "QC" , log_signal_stand:= (log_signal - log_signal_mean)/log_signal_sd ]
  
  # out
  fill.out$data <- dat
  return(fill.out)
  
}