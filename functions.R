require(readxl)
require(dplyr)
require(stringr)
require(tidyr)
require(data.table)
require(ggplot2)
require(MASS)
require(lme4)
require(car)
require(multcomp)
require(digest)
require(xlsx)

#### data loading and preprocessing #####
# Read the data from a MarkerLynx/EZInfo output .txt file
read.fun   <- function( file.path, meta.file, polarity) {
  
  #read raw data
  dat.raw  <- read.csv(file.path, sep = "\t", stringsAsFactors = F)
  
  # read meta data
  dat.meta <-  read_xlsx( meta.file)
  colnames(dat.meta) <- make.names(colnames(dat.meta))
  dat.meta <- subset(dat.meta, Polarity == polarity & is.na(Exclude) )
  
  # select all columns containing features
  dat.feat <- grep("[X,N,P][0-9]{1}", colnames(dat.raw))
  
  # Separate rtmz values and intensity values,
  # delete previous meta data and rename 1st column
  dat.rtmz  <- dat.raw[  1:2  , dat.feat]
  dat.clean <- dat.raw[-(1:2) , c(1,dat.feat)]
  colnames(dat.clean)[1] <- "File.Name"
  
  # bind values with new meta data
  dat.full <- left_join(dat.meta , dat.clean, by = "File.Name")
  
  # select all columns containing features
  dat.feat <- grep("[X,N,P][0-9]{1}", colnames(dat.full))
  
  # extract feature names, replace with accurate polarity identifier and update column names
  feat.names <- colnames(dat.full[ , dat.feat])
  feat.names <- str_replace(feat.names , "X", ifelse(polarity == "POS", "P", "N") )
  colnames(dat.full)[dat.feat] <- feat.names
  colnames(dat.rtmz) <- feat.names
  
  # transform to long format in data.table 
  dat.long <- data.table(gather(dat.full,key = "Feature", value = "Signal", dat.feat ))
  setkey(dat.long, Feature)
  
  # bind long data with rtmz data
  dat.rtmz <- data.table(Feature = colnames(dat.rtmz),
                         RT      = unlist(c(dat.rtmz[1,])),
                         mz      = unlist(c(dat.rtmz[2,])) )
  setkey(dat.rtmz, Feature)
  dat.long <- dat.long[dat.rtmz]
  
  return(dat.long)
}

# Clean data 
clean.fun  <- function( read.output, pars.file, polarity ) {
  
  # takes data from the read.fun() output
  
  ######################################################################################
  # This function applies several data cleaning steps to the raw auto-integrated data. #
  ######################################################################################
  
  # read parameter file with parameters for the different data cleaning steps
  pars <-  read_xlsx(pars.file, col_names = T, sheet = polarity )
  
  ######################### Heuristic denoise filter ########################################### 
  # The noise filter is based on the repetitions (lim_freq >5) of a                            #
  # *m/z* value within a certain mass window (lim_mz  > 0.02 Da) over a certain time-interval  #
  # (lim_rt >0.5 min) and not surpassing a certain intensity (im_sig < 10,000) The idea is     #
  # that if an *m/z* value repeats repeatedly over a long enough interval it is probably a     #
  # background ion coming from mobile phase.                                                   #
  ##############################################################################################
  
  # set thresholds for denoising
  lim_mz   <- pars$value[1]
  lim_freq <- pars$value[2]
  lim_rt   <- pars$value[3]
  lim_sig  <- pars$value[4]
  
  # Denoise data
  setkey(read.output, mz)
  dat.noise <- read.output[ , .(RT     = unique(RT),
                                mz     = unique(mz), 
                                Signal = unique(max(Signal, na.rm = T))) , by = Feature ]
  dat.noise[                                            , delta_mz    := c(0,diff(mz))                              ] 
  dat.noise[delta_mz < lim_mz                           , mz_idx1     := .I                                         ]
  dat.noise[delta_mz < lim_mz                           , mz_idx2     := .I                        , by = Feature   ]  
  dat.noise[                                            , noise_grp   := mz_idx2 - mz_idx1                          ] 
  dat.noise[!is.na(mz_idx1)                             , noise_freq  := .N                        , by = noise_grp ]
  dat.noise[noise_freq  >= lim_freq                     , delta_RT    := max(RT) - min(RT)         , by = noise_grp ]
  dat.noise[!is.na(noise_grp)                           , delta_sig   := max(Signal) - min(Signal) , by = noise_grp ] 
  dat.noise[delta_RT >= lim_rt & delta_sig  <= lim_sig  , noise       := 1L                        , by = noise_grp ]
  dat.noise <- dat.noise[noise == 1L]
  dat.noise[ ,max_sig := max(Signal),by = Feature]
  feat.filter_noise <- dat.noise$Feature
  
  # Plot noise features
  plot <- ggplot(dat.noise, aes(x = RT, y = mz) ) +
    geom_point(aes(col = log10(max_sig)), size = 1, alpha  = .8 ) +
    scale_color_gradient(low = "yellow",high = "red") +
    labs( title = "Deleted noise features",
          subtitle = paste0("delta mz < ", lim_mz,"; frequency > ", lim_freq, "; delta RT > ", lim_rt, ": max signal > ", lim_sig)) +
    theme_bw() +
    theme(legend.position = "bottom" )
  
  
  # Delete noise features
  read.output <- read.output[!(Feature %in% feat.filter_noise)]
  
  # Delete features with less than lim_qc (default = 2) values in qcs
  lim_qc <- pars$value[5]
  read.output[Signal == 0, Signal := NA]
  read.output[ , Delete := ifelse( sum(!is.na(Signal[Sample.Group == "QC"]) ) <= lim_qc, T, F ) , by = Feature ]
  feat.filter_qc <- unique(read.output[Delete == T, Feature])
  read.output <- read.output[!(Feature %in% feat.filter_qc)]
  
  # Delete features where each individual sample group, except QCs amd I_D_05, contain less than 50% of data available 
  lim_n_sig <- pars$value[6] 
  read.output[ Filler == "x" , Signal := 0 ] # fill [Signal] of filler-samples with 0 
  
  read.output[!(Sample.Group == "QC"), n_sig := sum(!is.na(Signal) ), by = c("Feature", "Sample.Group") ]
  read.output[ , Delete :=  { 
    t <- n_sig > lim_n_sig
    d <- ifelse(sum(t, na.rm = T) == 0, T, F )
    list(d)
  }, by = Feature ]
  feat.filter_grp<- unique(read.output[Delete == T, Feature])
  read.output <- read.output[!(Feature %in% feat.filter_grp)]
  
  # Delete features with more than "lim_na_glob" values missing over all sample groups except QCs/I_D_05
  lim_na_glob <- pars$value[7]
  read.output[ !(Sample.Group == "QC") , Delete := sum( is.na(Signal) ) >= lim_na_glob, by = Feature ]
  feat.filter_na <- unique(read.output[Delete == T, Feature])
  read.output <- read.output[!(Feature %in% feat.filter_na)]
  feat.included <-  unique(read.output$Feature)
  
  read.output[ Filler == "x" , Signal := NA ] # fill [Signal] of filler-samples with NA 
  
  return(list(data= read.output, 
              feature_data = list(
                deleted  = list(plot_noise   = plot,
                                noise_filter = feat.filter_noise,
                                filter_qc    = feat.filter_qc,
                                filter_grp   = feat.filter_grp,
                                filter_na    = feat.filter_na),
                included = feat.included) )
  )
}

# Median fold change adjustment
mfc.fun    <- function( clean.output ) {
  
  # takes data from the clean.fun() output
  
  ################################# median fold change (MFC) normalization #####################################
  # In order to make proper comparisons between samples, feature signals should be normalized between samples. #
  # Therefore we use the *median fold-change* method (MFC normalization). In short, a fold-change (FC) for     #
  # features per sample is calculated with respect to a reference sample. The refence sample is chosen to      #
  # be the one with the least missing values. Then, per sample, the median is taken over all the FC features.  #
  # This value is used to normalize all the features within a sample by deviding it by the Feature signal.     #    
  ##############################################################################################################
  
  # feature signal data
  dat <- clean.output$data
  
  # MFC calculations
  mfc.ref <- dat[Sample.Group != "QC", sum(!is.na(Signal) ), by = Sample.ID ][V1 == max(V1) , Sample.ID] # Select sample with least NA's, QCs excluded
  mfc.ref <- mfc.ref[1]                                                                                  # Choose 1st sample in case of ties
  dat[ , fold_change := Signal/Signal[Sample.ID == mfc.ref], by = Feature ]                              # Calulate fold change per feature
  dat[ ,  mfc_factor := median(fold_change, na.rm = T), by = Sample.ID]                                  # Calculte median fold change  for each sample
  dat[ ,  Signal_mfc := Signal/mfc_factor]                                                               # Devide signal by mfc factor
  
  # return mfc corrected signals and reference sample ID
  clean.output$data <- dat
  clean.output$mfc_ref <- mfc.ref
  
  return(clean.output)

}

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

#### merge, fill and standardize data ####  and log/standardize adjusted data ####
# merge positive and negative data
merge.fun <- function( dat.pos, dat.neg ) {
  
  dpos <- dat.pos$data
  dneg <- dat.neg$data
  
  # add ion mode
  dpos$mode <- "POS"
  dneg$mode <- "NEG"
  
  # mereg pos and neg
  dmerge <- rbind(dpos,dneg)
  
  # out
  out <- list(data = dmerge,
              features = unique(dmerge$Feature),
              meta = list(pos = dat.pos$feature_data, 
                          neg = dat.neg$feature_data))
  
  return (out)
  
}

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

#### model and feature selection ###
# model features with linear mixed effects model
mod.fun <- function( stnd.out, pars.file ) { #fill.output
  
  # data
  pars <-  read_xlsx(pars.file, col_names = T, sheet = "common" )
  dat <- stnd.out$data
  fts <- stnd.out$features
  
  # refactor
  grp <- unique(sort(dat$Sample.Group))[c(10,1:9,11:17)]                      # set levels, M_C_03 is reference group
  dat$Sample.Group <- factor(dat$Sample.Group, levels = grp)                  # relevel
  dat[Sample.Group != "QC" , SxC := paste0(Subject.Species, "_", Challange) ] # make new groups Species x Challenge
  
  # models
  n <- length(fts)
  mods <- dat[Sample.Group != "QC", {
    f <- unique(Feature)
    cat(f," ", which(fts == f), " from ", n , "features  ") #counter
    m <- list(lmer(log_signal_stand ~ Sample.Group + (1  | SxC : Animal), REML = T ) ) #model
    cat("\n")
    list(m)
  }, by = Feature]
  setnames(x = mods, old = "m", new = "lmer_mods")
  
  # get p vals
  cat("\nCalculate P values.\n")
  alpha_lmer <- pars$value[11]
  mods[  , p_lmer := sapply(lmer_mods, function(m) Anova(m)["Pr(>Chisq)"] ), by = Feature ]
  
  # check if model is singular ( see ?isSingular() )
  mods[  , singular := sapply(lmer_mods, function(m) isSingular(m) ), by = Feature ]
  
  # out
  stnd.out$models <- mods
  return(stnd.out)
  
}

# Post-hoc analysis on relevant contrasts
posthoc.fun <- function( mod.out, contrast.file ) { #mod.output
  
  mods <- mod.out$models
  K    <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  fts    <- mod.out$models[p_lmer < 5e-3, Feature] 
  n    <- length(fts)
  
  mods.ph <- sapply( fts , function(f) {
    cat(f," ", which(fts == f), " from ", n , "features\n") #counter
    sapply(mods[Feature == f]$lmer_mods, function(m) {
      summary(glht(m , linfct = K ) )
      }, simplify = F )
    }, simplify = F, USE.NAMES = T )
  
  #out
  mod.out$posthoc <- mods.ph
  return(mod.out)
  
}

# Select features
select.fun <- function( posthoc.out, pars.file, contrast.file  ) { #posthoc.output
  
  # load parameter data
  pars <-  read_xlsx(pars.file, col_names = T, sheet = "common" )
  
  # set params
  alpha_select <- 0.05 #pars$value[12]
  mods <- posthoc.out$posthoc
  f    <- posthoc.out$models[p_lmer < 5e-3, Feature] #features[1:3]  #
  K    <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  
  # extract pvals from posthoc results
  m_extr <- sapply(f, function(f) 
    sapply(mods[[f]], function(m1) unlist(m1, recursive = F),
           simplify = F, USE.NAMES = T) , simplify = F, USE.NAMES = T )
  
  p_vals      <- data.frame(sapply(f, function(f) m_extr[[f]][[1]]$test.pvalues ))
  p_vals$grps <- rownames(K)
  p_vals      <- data.table(gather(p_vals,key = Feature, value = pvals, 1:length(f)))
  
  # select C57 M vs I
  p_vals[  , C57_MI :=  ifelse(sum(pvals[1] < alpha_select | pvals[2] < alpha_select) == 1 & sum(pvals[3:6] < alpha_select) >= 1, T,  F) , by = Feature]
  f_C57 <- p_vals[C57_MI == T, unique(Feature)]
  
  # select DBA M vs I
  p_vals[  , DBA_MI  :=  ifelse( sum( pvals[7:8] <  alpha_select) >= 1, T,  F) , by = Feature]
  f_DBA <- p_vals[DBA_MI == T, unique(Feature)]
  
  # select DBA vs C57
  p_vals[  , C57_DBA :=  ifelse( sum( pvals[9:11] < alpha_select) >= 2, T,  F) , by = Feature]
  f_IS <- p_vals[C57_DBA == T, unique(Feature)]
  
  f_total <- unique(c(f_C57,f_DBA,f_IS))
  
  # out
  out <- list( data = posthoc.out$data,
               features = posthoc.out$features,
               models = posthoc.out$models,
               posthoc = posthoc.out$posthoc,
               meta = posthoc.out$meta,
               selections = list( pvals = p_vals,
                                  select_C57 = f_C57,
                                  select_DBA = f_DBA,
                                  select_InterSpecies = f_IS,
                                  select_total = f_total))
  
  return(out)
  
}



