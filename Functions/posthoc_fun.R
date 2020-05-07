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