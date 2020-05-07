# bind manually reintegrated TargetLynx data to auto-itegrated data
reint.bind.fun <- function( select.output, meta.file, reint.files ) {
  
  # functions
  feature_extract <- function(d) {
    f <- str_extract(d$V1,"[NPX]\\d{1,}\\.\\d{1,}\\_\\d{2,}\\.\\d{3,}.*")
    d <- data.table(Feature = f, Polarity = unique(d$Polarity))
    d <- d[!is.na(Feature)]
    d[ , Feature := ifelse(Polarity == "POS", str_replace(Feature,"^X","P"),str_replace(Feature,"^X","N"))]
    return(d$Feature)
  }
  
  # read meta data, make propper col names, exclude samples
  dat.meta <-  read_xlsx(meta.file)
  col.names.meta <- make.names(colnames(dat.meta) )
  colnames(dat.meta) <- col.names.meta
  dat.meta <- subset(dat.meta, is.na(Exclude) )
  dat.meta <- data.table(dat.meta)
  col.names.meta <- colnames(dat.meta)
  
  # Read TargetLynx results.txt files
  dat <- sapply(reint.files, function(file) {
    d <- read.csv(file , sep = "\t", header = FALSE, stringsAsFactors = F) 
    p <- ifelse(grepl("POS_REINT",file),"POS","NEG") # polarity set by file name!
    d <- data.table(d, Polarity = p)
  }, USE.NAMES = T,simplify = F)
  
  # extract metabolite names from all reintegrations and set polarity indicator (N/P)
  f_all <- sapply(dat, function(d) feature_extract(d), USE.NAMES = F, simplify = T)
  f_all <- do.call(unlist, list(f_all, use.names = F)) 
  
  # check if all selected raw features are present in reintegration data, if not --> reintegrate the ones that are missing
  f_sel  <- select.output$selections$select_total
  
  if ( all(f_sel %in% f_all) ) {
    
    # read and wrangle TargetLynx data
    dat <- sapply(dat, function(d) {
      
      # extract metabolite names from individual QL files and set propper polarity indicator (N/P)
      f_ind <- feature_extract(d)
      
      # determine if .txt file was resaved in Excel (causes blank rows to be removed and changes row numbers)
      # and determine sample numbers, column names and sample rows
      resaved <- ifelse(d[2,1] == "", T, F)
      if (resaved == T){
        # #For "re-saved TargetLynx text files
        ql.id     <- d$V1[3]
        n.samples <- max(as.integer(d$V2[grep("[0-9]",d$V2)]))
        c.names   <- make.names(as.character(d[7, -(1:2) ]))
        r.select  <- as.vector(sapply(seq_along(f_ind), function(i) (8 + (i - 1) * (n.samples + 4) ) : (3 + i * (n.samples + 4)  ) ) )
      } else {
        #For untouched TargetLynx text files
        ql.id     <- d$V1[2]
        n.samples <- max(as.integer(d$V2[grep("[0-9]",d$V2)]))
        c.names   <- make.names(as.character(d[4, -(1:2) ]))
        r.select  <- as.vector(sapply(seq_along(f_ind), function(i) (5 + (i - 1) * (n.samples + 2) ) : (2 + i * (n.samples + 2)  ) ) )
      }
      
      # adjsut colnames 
      c.names <- ifelse(c.names %in%  c("POS","NEG"), "Polarity", c.names)
      
      # tidy up TargetLynx data
      d <- d[r.select, -(1:2) ]
      colnames(d)    <- c.names
      colnames(d)[1] <- "File.Name"
      d$Area               <- as.numeric(d$Area)
      d$RT                 <- as.numeric(d$RT)
      d$Peak.Start.Height  <- as.numeric(d$Peak.Start.Height)
      d$Height             <- as.numeric(d$Height)
      d$S.N                <- as.numeric(d$S.N)
      d$Feature            <- as.vector(sapply(f_ind, function(m)  rep(m,n.samples) ) )
      d$QL_ID              <- ql.id
      
      d <- data.table(d)
      setnames(d,"RT","RT_reint_QL")
      dat.meta <- data.table(dat.meta)
      
      pol <- unique(d$Polarity)
      dat.meta <- dat.meta[Polarity == pol]
      d_rt <- d[ , .(Feature, File.Name, RT_reint_QL)]
      
      d <- dcast(d, File.Name + Sample.Text + Vial + QL_ID  ~ Feature, value.var = "Area")
      d <- merge(dat.meta, d, by = "File.Name", all.x = T, variable.factor = F)  #  # check!!! --> this should automatically delete excluded samples
      d[ , QL_ID := unique(na.omit(QL_ID))] # to give filler samples correct QL_ID in order to remove duplicates later
      d <- melt(d,measure.vars = f_ind)
      setnames(d,c("variable", "value"),c("Feature", "Area_reint"))
      
      d <- merge(d, d_rt, by = c("File.Name", "Feature"), all.x = T, variable.factor = F)
      
      return(d)
      
    }, simplify = F )
    dat <- do.call(rbind, dat)
    
    # delete superfluous features
    dat[ , idx := as.integer(factor(QL_ID)), by = Feature] # check for duplicates
    dat <- dat[idx == 1L,]                                 # delete duplicates
    dat <- dat[!grepl("[*]",Feature)]                      # delete features tagged with an asterix
    
    # cross-check reintegrated features with cleaned features and delete if necesary
    f_deleted <- do.call(unlist,list(select.output$meta$pos$deleted, 
                                     select.output$meta$neg$deleted, 
                                     use.names = F))            # features that were deleted from raw data
    f_reint   <- unique(dat$Feature)                            # all reintegrated features
    f_dat     <- f_reint[!(f_reint %in% f_deleted)]             # get remaining features for melting dat later
    dat <- dat[Feature %in% f_dat] 
    
    # calculate adjusted RT by taking the weighted average of RTs over samples of reintegrated features
    dat[ , RT_reint := weighted.mean(RT_reint_QL,Area_reint, na.rm = T), by = Feature]
    
    # try to extract RTMZ values from Feature names otherwise assign NA
    mzrt <- str_extract_all(f_dat, "[0-9]{1,4}[.][0-9]{2,4}" )
    mzrt <- data.table(t(sapply(mzrt, function(m) {
      m <- as.double(m)
      if (length(m) == 0 ) { m <- c(NA,NA) } else { m <- m }
      return(m)
    }, simplify = T, USE.NAMES = T )))
    setnames(mzrt, c("V1", "V2"), c("RT_extr", "mz_reint" ) )
    mzrt$Feature <- f_dat
    dat <- merge(dat, mzrt, by = "Feature" , variable.factor = F)
    
    # set reintegrated flag
    dat$reintegrated <- T
    
    # merge raw data with reintegrated data  
    dat <- merge(select.output$data , dat, by = c( "Feature", col.names.meta), all = T , variable.factor = F)
    dat <- dat[is.na(reintegrated), reintegrated := F]
    
    # prepare output
    select.output$data <- dat
    select.output$feature_data$reint_all <- f_dat
    
    return(select.output)
    
  } else {
    
    # print warnings and features missing from the reintegrations
    print("Not all selected features are reintegrated. Reintegrate the following:")
    print(f_sel[!(f_sel %in% f_all)])
    
    return(select.output)
    
  }
}