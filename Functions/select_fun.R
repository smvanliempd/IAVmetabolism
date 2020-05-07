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
  
  # select C57 Mock vs Infected
  p_vals[  , C57_MI :=  ifelse(sum(pvals[1] < alpha_select | pvals[2] < alpha_select) == 1 & sum(pvals[3:6] < alpha_select) >= 1, T,  F) , by = Feature]
  f_C57 <- p_vals[C57_MI == T, unique(Feature)]
  
  # select DBA Mock vs Infected
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