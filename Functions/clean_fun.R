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
