# Select features
select.fun <- function( posthoc.out, pars.file, contrast.file  ) {
  
  # load data
  pars <- read_xlsx(pars.file, col_names = T, sheet = "COMMON" )
  K    <- as.matrix(read.csv( contrast.file, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE, row.names = 1) )
  mods <- posthoc.out$posthoc
  
  # set feature selection params (alphas)
  # alpha_lmer   <- pars$value[11]
  alpha_select <- pars$value[2] #12
  
  # get features that were selected for post-hoc analysis
  fts    <- names(dat.mod$posthoc) #posthoc.out$models[p_lmer < alpha_lmer, Feature]
  
  # extract pvals from posthoc results
  m_extr <- sapply(fts, function(f) 
    sapply(mods[[f]], function(m1) unlist(m1, recursive = F),
           simplify = F, USE.NAMES = T) , simplify = F, USE.NAMES = T )
  p_vals      <- data.frame(sapply(fts, function(f) m_extr[[f]][[1]]$test.pvalues ))
  p_vals$grps <- rownames(K)
  p_vals      <- data.table(gather(p_vals,key = Feature, value = pvals, 1:length(fts)))
  
  # select C57 Mock vs Infected
  p_vals[  , C57_MI :=  ifelse(sum(pvals[1] < alpha_select | pvals[2] < alpha_select) == 1 & sum(pvals[3:6] < alpha_select) >= 1, T,  F) , by = Feature]
  f_C57 <- p_vals[C57_MI == T, unique(Feature)]
  
  # select DBA Mock vs Infected
  p_vals[  , DBA_MI  :=  ifelse( sum( pvals[7:8] <  alpha_select) >= 1, T,  F) , by = Feature]
  f_DBA <- p_vals[DBA_MI == T, unique(Feature)]
  
  # select DBA vs C57
  p_vals[  , C57_DBA :=  ifelse( sum( pvals[9:11] < alpha_select) >= 2, T,  F) , by = Feature]
  f_IS <- p_vals[C57_DBA == T, unique(Feature)]
  
  # totoal amount of markers
  f_total <- unique(c(f_C57,f_DBA,f_IS))
  
  # out
  posthoc.out$selections <- list( pvals = p_vals,
                                 select_C57 = f_C57,
                                 select_DBA = f_DBA,
                                 select_InterSpecies = f_IS,
                                 select_total = f_total)

  
  # out <- list( data = posthoc.out$data,
  #              features = posthoc.out$features,
  #              models = posthoc.out$models,
  #              posthoc = posthoc.out$posthoc,
  #              meta = posthoc.out$meta,
  #              selections = list( pvals = p_vals,
  #                                 select_C57 = f_C57,
  #                                 select_DBA = f_DBA,
  #                                 select_InterSpecies = f_IS,
  #                                 select_total = f_total))
  
  return(posthoc.out)
  
}
