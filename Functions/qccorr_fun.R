# Quality control correction
qccorr.fun <- function( mfc.output, pars.file, polarity ) {
  
  # Takes output from mfc.fun()
  
  ######################### Quality control corrections #############################
  # If features show drift in quality control (QC) samples after MFC normalization, #
  # then this drift will be corrected with a polynomial least square fitted model   #
  # with the injection number as the idependent variable and the scaled QC signal   #
  # as the dependend variable.                                                      #
  ###################################################################################
  
  # read parameter file
  pars <-  read_xlsx(pars.file, col_names = T, sheet = polarity )
  
  # input data and parameters
  dat          <- mfc.output$data
  lim_qccorr   <- pars$value[8]
  alpha_qccorr <- pars$value[9]
  pol.order    <- pars$value[13]
  
  # Data quality check and transformations
  dat[Sample.Group == "QC" , qc_corr  := sum(!is.na(Signal_mfc) ) >= lim_qccorr , by = Feature ] # determine if there are enough QC samples
  dat[Sample.Group == "QC" & !is.na(Signal_mfc) , qc_ref_signal := Signal_mfc[1] , by = Feature ]    # determine first non-NA QC signal -> reference value
  dat[qc_corr == T         , sig_scld := Signal_mfc/qc_ref_signal  , by = Feature ]              # make scaled QC signals based on reference value
  dat[   , c("qc_corr_fact","qc_corr_flag") := {     
    
    # Get QC data
    d <- data.frame( inj = Injection.Number,  
                     sig = sig_scld)
    d <- d[complete.cases(d), ]
    
    # If there are more than lim_qccorr data point in QC samples make QC samples
    if (nrow(d) >= lim_qccorr )  { 
      
      # model with parameterized polynomial order 
      # (should probably not be higher than 3 (cubic) ).
      # Models are accepted if p_corr < alpha_qccorr for at least 1 coefficient.
      o     <- 1:pol.order
      fun   <- formula(paste0("sig ~ ",paste(paste0("I(inj^", o, ")" ), collapse  = "+") ) )
      mod   <- summary(rlm(fun, data = d))
      
      # model summary
      tvals  <- mod$coefficients[,"t value"]
      df     <- mod$df[2] 
      pvals  <- dt(tvals, df)
      
      # at least one of the pvals should be < alpha except for intercept
      if (sum(pvals[2:(pol.order+1)] < alpha_qccorr) > 0 ) { 
        
        cf <- mod$coefficients[,"Value"]
        qc_fact <- sapply( 1:(pol.order+1) , function(i) cf[i] * Injection.Number^(i - 1), simplify = T  )
        qc_fact <- rowSums(qc_fact)
        qc_flag <- T
        
      } else {
        
        qc_fact <- 1
        qc_flag <- F
        
      }
    } else {
      
      qc_fact <- 1
      qc_flag <- F
      
    }
    
    list(qc_fact, qc_flag )
    
  }, by = Feature ]
  
  feat.qccorr <- unique(dat[qc_corr_flag == T, Feature ])
  dat[qc_corr_fact < 0.1, qc_corr_fact := 0.1 ] # all corrections <10-fold set to 10-fold
  dat[ , Signal_mfc_qc := ifelse(qc_corr_flag == T, Signal_mfc/qc_corr_fact, Signal_mfc)]
  
  # prepare output data
  mfc.output$data <- dat
  mfc.output$feature_data$qc_corr <- feat.qccorr
  return(mfc.output)
  
}