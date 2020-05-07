# Clean reintegrated data
reint.clean.fun <- function( reint.output, pars.file ) {
  
  # read data
  dat.clean <- reint.output$data
  f_reint   <- reint.output$feature_data$reint_all
  
  for(plr in c("POS","NEG") ) {
    
    # read parameter file
    pars <-  read_xlsx(pars.file, col_names = T, sheet = plr )
    
    # Delete features with less than 2 values in qcs
    lim_qc <- pars$value[5]
    dat.clean[reintegrated == T & Polarity == plr, Delete := ifelse( sum(!is.na(Area_reint[Sample.Group == "QC"]) ) <= lim_qc, T, F ) , by = Feature ]
    feat.filter_qc <- unique(dat.clean[Delete == T & Polarity == plr, Feature])
    dat.clean <- dat.clean[!(Feature %in% feat.filter_qc)]
    
    # Delete features where each individual sample group, except QCs amd I_D_05, contain less than 50% of data available
    lim_n_sig <- pars$value[6]
    dat.clean[ Filler == "x" & reintegrated == T & Polarity == plr, Area_reint := 0 ] # fill [Area_reint] of filler-samples with 0
    dat.clean[!(Sample.Group == "QC") & reintegrated == T & Polarity == plr, n_sig := sum(!is.na(Area_reint) ), by = c("Feature", "Sample.Group") ]
    dat.clean[ reintegrated == T & Polarity == plr, Delete :=  {
      t <- n_sig > lim_n_sig
      d <- ifelse(sum(t, na.rm = T) == 0, T, F )
      list(d)
    }, by = Feature ]
    feat.filter_grp <- unique(dat.clean[Delete == T & Polarity == plr, Feature])
    dat.clean <- dat.clean[!(Feature %in% feat.filter_grp)]
    
    # Delete features with more than "lim_na_glob" values missing over all sample groups except QCs/I_D_05
    lim_na_glob <- pars$value[7]
    dat.clean[ !(Sample.Group == "QC") & reintegrated == T & Polarity == plr, Delete := sum( is.na(Area_reint) ) >= lim_na_glob, by = Feature ]
    feat.filter_na <- unique(dat.clean[Delete == T & Polarity == plr, Feature])
    dat.clean <- dat.clean[!(Feature %in% feat.filter_na)]
    
    if (plr == "POS") {
      reint.output$meta$pos$deleted$filter_qc_reint   <- feat.filter_qc
      reint.output$meta$pos$deleted$filter_grp_reint <- feat.filter_grp
      reint.output$meta$pos$deleted$filter_na_reint  <- feat.filter_na
    } else {
      reint.output$meta$neg$deleted$filter_qc_reint   <- feat.filter_qc
      reint.output$meta$neg$deleted$filter_grp_reint <- feat.filter_grp
      reint.output$meta$neg$deleted$filter_na_reint  <- feat.filter_na
    }
    
  }
  
  # get remaining features
  feat.included <-  unique(dat.clean[reintegrated == T ,Feature])
  reint.output$feature_data$reint_incl <- feat.included
  
  # fill [Area_reint] of filler-samples with NA
  dat.clean[ Filler == "x" , Area_reint := NA ]
  
  # out
  reint.output$data <- dat.clean
  return(reint.output)
  
}