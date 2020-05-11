# select reintegrated feautures
reint.select.fun <- function( reint_ph.output, pars.file, contrast.file  ) {
  
  # load parameter data
  pars <-  read_xlsx(pars.file, col_names = T, sheet = "COMMON" )
  alpha_select <- pars$value[14]
  
  # load data
  mods <- reint_ph.output$posthoc
  f    <- reint_ph.output$feature_data$reint_incl # select all included reintegrated features
  K    <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  
  # extract pvals from posthoc results
  m_extr <- sapply(f, function(f) 
    sapply(mods[[f]], function(m) unlist(m, recursive = F), 
           simplify = F, USE.NAMES = T) , simplify = F, USE.NAMES = T )
  
  
  p_vals      <- data.table(sapply(f, function(f) m_extr[[f]][[1]]$test.pvalues ))
  p_vals$grps <- rownames(K)
  p_vals <- melt(p_vals,  measure.vars = f, variable.name = "Feature", value.name = "pvals", variable.factor = F )
  
  # select C57 M vs I
  p_vals[  , C57_MI :=  ifelse( sum( pvals[1:6] <  alpha_select ) >= 1, T,  F) , by = Feature]
  f_C57 <- p_vals[C57_MI == T, unique(Feature)]
  
  # select DBA M vs I
  p_vals[  , DBA_MI  :=  ifelse( sum( pvals[7:8] <  alpha_select) >= 1, T,  F) , by = Feature]
  f_DBA <- p_vals[DBA_MI == T, unique(Feature)]
  
  # select DBA vs C57
  p_vals[  , C57_DBA :=  ifelse( sum( pvals[9:11] < alpha_select) >= 2, T,  F) , by = Feature]
  f_IS <- p_vals[C57_DBA == T, unique(Feature)]
  
  feat_total <- unique(c(f_C57,f_DBA,f_IS))
  
  reint_ph.output$selections_reint <- list(pvals = p_vals,
                                           select_C57 = f_C57,
                                           select_DBA = f_DBA,
                                           select_InterSpecies = f_IS,
                                           select_total = feat_total)
  
  return(reint_ph.output)
  
}
