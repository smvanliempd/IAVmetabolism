# model features with linear mixed effects model
mod.fun <- function( stnd.out ) {
  
  # takes log/scaled output
  
  # data
  dat <- stnd.out$data
  fts <- stnd.out$features
  
  # refactor
  grp <- unique(sort(dat$Sample.Group))[c(10,1:9,11:17)]                      # set levels, M_C_03 is reference group
  dat$Sample.Group <- factor(dat$Sample.Group, levels = grp)                  # relevel
  dat[Sample.Group != "QC" , SxC := paste0(Subject.Species, "_", Challenge) ] # make new groups Species x Challenge
  
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
  mods[  , p_lmer := sapply(lmer_mods, function(m) Anova(m)["Pr(>Chisq)"] ), by = Feature ]
  
  # check if model is singular ( see ?isSingular() )
  mods[  , singular := sapply(lmer_mods, function(m) isSingular(m) ), by = Feature ]
  
  # out
  stnd.out$models <- mods
  return(stnd.out)
  
}