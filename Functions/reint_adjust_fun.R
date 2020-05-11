# MFC/QC adjust reintegrated data
reint.adjust.fun <- function( reint.clean.output, pars.file) {
  
  # read data and vars
  dat          <- reint.clean.output$data
  f_reint      <- reint.clean.output$feature_data$reint_incl
  
  ## MFC adjustments
  # Fill MFC values for reintegrations
  dat[ ,mfc_factor := mean(mfc_factor, na.rm = T), by = File.Name] 
  
  # make MFC adjustment
  dat[Feature %in% f_reint , Area_reint_mfc := Area_reint/mfc_factor]
  
  ## QC adjustments for reintegrations in positive and negative ionization
  for(plr in c("POS","NEG") ) {
    
    # read parameter file
    pars <-  read_xlsx(pars.file, col_names = T, sheet = plr )
    lim_qccorr   <- pars$value[8]
    alpha_qccorr <- pars$value[9]
    pol.order    <- pars$value[13]
    
    ## data quality check and scaling
    # determine if there are enough QC samples
    dat[Sample.Group == "QC" & Feature %in% f_reint & Polarity == plr, qc_corr_reint  := sum(!is.na(Area_reint_mfc) ) >= lim_qccorr , by = Feature ]
    
    # determine first non-NA QC signal -> reference value
    dat[Sample.Group == "QC" & !is.na(Area_reint_mfc) & Feature %in% f_reint & Polarity == plr, qc_ref_Area_reint := Area_reint_mfc[1] , by = Feature ]   
    
    # make scaled QC signals based on reference value
    dat[qc_corr_reint == T  & Feature %in% f_reint & Polarity == plr, sig_scld_reint := Area_reint_mfc/qc_ref_Area_reint  , by = Feature ]              
    
    # make qc models
    dat[Feature %in% f_reint & Polarity == plr, c("qc_corr_fact_reint","qc_corr_flag_reint") := {     
      
      #data
      d <- data.frame( inj = Injection.Number, sig = sig_scld_reint)
      d <- d[complete.cases(d), ]
      
      if (nrow(d) >= lim_qccorr )  { # check if there are enough qc samples
        
        #model
        o     <- 1:pol.order
        fun   <- formula(paste0("sig ~ ",paste(paste0("I(inj^", o, ")" ), collapse  = "+") ) )
        mod   <- summary(rlm(fun, data = d))
        
        # summary
        tvals  <- mod$coefficients[,"t value"]
        df     <- mod$df[2] 
        pvals  <- dt(tvals, df)
        
        if (sum(pvals[2:(pol.order+1)] < alpha_qccorr) > 0 ) { # at least one of the pvals should be < alpha except for intercept
          
          c <- mod$coefficients[,"Value"]
          qc_fact <- sapply( 1:(pol.order+1) , function(i) c[i] * Injection.Number^(i - 1), simplify = T  )
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
    
    # qc correction
    dat[qc_corr_fact_reint < 0.1 & Polarity == plr, qc_corr_fact_reint := 0.1 ] # all corrections <10-fold set to 10-fold
    dat[Polarity == plr, Area_reint_mfc_qc := ifelse(qc_corr_flag_reint == T, Area_reint_mfc/qc_corr_fact_reint, Area_reint_mfc)]
    
    # which reintegrated features were qc corrected?
    feat.qccorr <- unique(dat[qc_corr_flag_reint == T & Polarity == plr, Feature ])
    
    # Collect QC corrected features 
    if (plr == "POS") {
      reint.clean.output$meta$pos$qc_corr_reint <- feat.qccorr 
    } else {
      reint.clean.output$meta$neg$qc_corr_reint <- feat.qccorr 
    }
    
  }
  
  # out
  reint.clean.output$data <- dat
  return(reint.clean.output)
  
}