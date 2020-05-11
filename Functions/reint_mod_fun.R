# model and select reintegrated data
reint.mod.fun <- function( stnd.out, pars.file, contrast.file ) {
  
  # data
  K <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  pars <-  read_xlsx(pars.file, col_names = T, sheet = "COMMON" )
  dat <- stnd.out$data
  f_incl <- stnd.out$feature_data$reint_incl
  n <- length(f_incl)
  
  # refactor
  grp <- unique(sort(dat$Sample.Group))[c(10,1:9,11:17)]                      # set levels, M_C_03 is reference group
  dat$Sample.Group <- factor(dat$Sample.Group, levels = grp)                  # relevel
  dat[Sample.Group != "QC" , SxC := paste0(Subject.Species, "_", Challenge) ] # make new groups Species x Challenge
  
  # models
  cat("Make LME models./n/n")
  mods <- dat[Feature %in% f_incl & Sample.Group != "QC", {
    f <- unique(Feature)
    cat("LME model for ",f," : ", which(f_incl == f), " from ", n , "features  ") #counter
    m <- list(lmer(reint_log_stand ~ Sample.Group + (1  | SxC : Animal), REML = T ) )
    cat("\n")
    list(m)
  }, by = Feature]
  setnames(x = mods, old = "m", new = "lmer_mods")
  
  mods$V1
  
  # get p vals and select relevant models
  alpha_lmer <- pars$value[11] 
  mods[  , p_lmer := sapply(lmer_mods, function(m) Anova(m)["Pr(>Chisq)"] ), by = Feature ]
  
  # post hoc analysis
  cat("\nPost-hoc analysis.\n\n")
  mods.ph <- sapply( f_incl , function(f){
    cat("Post-hoc test for ",f," : ", which(f_incl == f), " from ", n , "features\n") #counter
    sapply(mods[Feature == f]$lmer_mods, function(m) {
      summary(glht(m , linfct = K ) )
    }, simplify = F )
  }, simplify = F, USE.NAMES = T )
  
  # out
  stnd.out$data <- dat
  stnd.out$models <- mods
  stnd.out$posthoc <- mods.ph
  
  return(stnd.out)
  
}
