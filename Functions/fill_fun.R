# fill empty sample groups with minimal signals (25% percentile) values with a 15% random error
fill.fun <- function( qccorr.output ) {
  
  dat <- qccorr.output$data
  
  qnt25       <- dat[ , quantile(.SD[ , min(Signal, na.rm = T), by = Feature ]$V1, probs = 0.25), ] # set minimum fill signal - 25% qunatile
  incomp.grp  <- unique(dat[Filler == "x", Sample.Group])                                           # determine which sample groups contain filler-samples
  fill.sample <- unique(dat[(Sample.Group %in% incomp.grp) & (is.na(Filler) ), Sample.ID])          # deteremine real samples in groups that contain filler-samples
  dat[ n_sig == 0  | (Sample.ID %in% fill.sample & is.na(Signal_mfc_qc) ),  fill := T]              # set fill flag 
  
  
  set.seed(19743009)
  dat[fill == T , Signal_mfc_qc := runif(.N, 0.85, 1.15) * qnt25  ] # fill flagged samples with [qnt25] ± 15%
  
  qccorr.output$data <- dat
  return(qccorr.output)
  
}