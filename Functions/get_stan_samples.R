get.stan.samples <- function(dat ,mod , n_chains = 4, get_delta_bl_samples = F) {
  
  # function to transform log2 data in a %difference bounded by -100%
  log2delta <- function(d) {  100 * ( (2^d)  -  1 ) }
  
  # seed for reproducibility
  seedy <- 20150613
  
  # stan parameters
  n_iter <- 5500
  n_warm <- 1500
  
  # features
  fts <- names(dat)
  n_fts <- length(fts)
  
  # quantiles for samples
  qtls <- c(5,25,50,75,95)/100
  
  # samples
  ss <- sapply(fts, function(f) {
    
    # Timing and mesages
    ptm <- proc.time()
    cat(f,"... ",which(fts == f)," of ",n_fts)
    
    # Samples to extract from the model
    sample.pars <- c("del_c_3i_0i",
                     "del_c_5i_0i",
                     "del_c_8i_0i",
                     "del_c_18i_0i",
                     "del_c_30i_0i",
                     "del_c_3i_3m",
                     "del_c_5i_5m",
                     "del_c_8i_8m",
                     "del_c_18i_18m",
                     "del_c_30i_30m",
                     "del_d_3i_0i",
                     "del_d_5i_0i",
                     "del_d_3i_3m",
                     "del_d_5i_5m",
                     "del_d3m_c3m",
                     "del_d5m_c5m",
                     "del_d0i_c0i",
                     "del_d3i_c3i",
                     "del_d5i_c5i",
                     "area_gr_hat" )
    
    # Samples to pool for inter-strain baseline differences
    delta.bl.pars <- c("del_d3m_c3m",
                       "del_d5m_c5m",
                       "del_d0i_c0i")
    
    # extract delta-samples
    s <- sampling(object = mod,
                  pars =  sample.pars,
                  data = dat[[f]],
                  warmup = n_warm ,
                  iter = n_iter ,
                  chains= n_chains ,
                  cores = n_chains,
                  control = list(adapt_delta = 0.99), 
                  seed = seedy )
    s <- extract(s,pars = sample.pars)
    s <- do.call(cbind,s)
    
    # change log2 values of differences to delta percentages
    # s[ , sample.pars[-20]] <- log2delta( s[ , sample.pars[-20]])
    
    # all samples for interstrain differences for baselines (t0, d3m, d5m)
    s_bl <- c(s[,delta.bl.pars])
    
    # sample 6000 times from inter-strain differences for baselines (t0, d3m, d5m)
    set.seed(seedy)
    s_pl_bl <- data.table(s_pl_bl = sample(x = s_bl, 
                                           size = 6000,
                                           replace = F))
    s_pl_bl$Feature <- f
    
    # calculate del_del_cd_3i0i, i.e. the difference in percent-change between strains
    #  these are *no* log2 values!!!
    del_del_str_inf <-  c(log2delta(s[,"del_d_3i_0i"])) -  log2delta(c(s[,"del_c_3i_0i"]))
    s <- cbind(s, del_del_str_inf)
    
    
    # calculate proportion of samples bigger than 0
    p <- data.table(t(apply(s,2, function(a) sum(a > 0 )/length(a))))
    p$Feature <- f
    p$metric  <- "P_pos"
    
    # calculate proportion of samples bigger than 0 for pooled baselines
    p_pl_bl <- sum(s_bl > 0 )/length(s_bl)
    
    # calulate quantiles
    q <- data.table(apply(s,2, function(a) quantile(a, qtls)))
    q$Feature <- f
    q$metric <- paste0("q",round(qtls*100,0))
    
    # calulate quantiles for pooled baselines
    q_pl_bl <- quantile(s_bl, qtls)
    
    # bin p and q
    d <- rbind(q,p)
    
    # add the pooled baseline p and q values to the rest of the data
    d$del_pl_bl <- c(q_pl_bl,p_pl_bl)
    
    # progress message
    cat( " done in ",c(proc.time() - ptm)[3]," sec.\n")
    
    # return
    return(list(d = d,s_pl_bl = s_pl_bl))
    
  }, USE.NAMES = T, simplify = F)
  
  s_quantiles <- sapply(ss, function(l) l$d, simplify = F,USE.NAMES = T)
  s_quantiles <- do.call(rbind, s_quantiles)
  s_quantiles <- melt(s_quantiles, id.vars = c("Feature","metric"),variable.factor = F)
  s_quantiles <- dcast(s_quantiles, ... ~ metric)
  
  # just in case you want to get a density plot of the log2 delta baseline samples
  if(get_delta_bl_samples == T) {
    
    delta_bl_samples <- sapply(ss, function(l) l$s_pl_bl, simplify = F,USE.NAMES = T)
    delta_bl_samples <- do.call(rbind, delta_bl_quantiles)
    return(list(s_quantiles=s_quantiles,delta_bl_samples=delta_bl_samples))
    
  } else {
    
    return(list(s_quantiles=s_quantiles)) 
    
  }
}
