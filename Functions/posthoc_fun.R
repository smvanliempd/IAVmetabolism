# Post-hoc analysis on relevant contrasts
posthoc.fun <- function( mod.out, pars.file, contrast.file ) {
  
  # load data
  pars <- read_xlsx(pars.file, col_names = T, sheet = "COMMON" )
  mods <- mod.out$models
  K    <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  
  # select lmer mods for posthoc tests based on alpha_lmer
  alpha_lmer <- pars$value[1] #11
  fts    <- mod.out$models[p_lmer < alpha_lmer, Feature] 
  n    <- length(fts)
  
  # posthoc tests
  mods.ph <- sapply( fts , function(f) {
    cat(f," ", which(fts == f), " from ", n , "features\n") #counter
    sapply(mods[Feature == f]$lmer_mods, function(m) {
      summary(glht(m , linfct = K ) )
    }, simplify = F )
  }, simplify = F, USE.NAMES = T )
  
  # out
  mod.out$posthoc <- mods.ph
  return(mod.out)
  
}