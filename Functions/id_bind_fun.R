# bind results for reintegrated features with previous identifications and export the feature-list to excel file.
id.bind.fun <- function( metb.cor.out, id.file ) {
  
  # input data
  dat    <- metb.cor.out$data
  f_C57  <- metb.cor.out$selections_reint$select_C57
  f_DBA  <- metb.cor.out$selections_reint$select_DBA
  f_IS   <- metb.cor.out$selections_reint$select_InterSpecies
  clstr  <- metb.cor.out$correlations$metabolic$`C57BL/6`$baseline$nodes[ , .(id, metb_cor_grp)]
  
  # Take ALL reintegrated features 
  f_rep    <- dat[ reintegrated == T , .(rep_feature = rep_feature[1], 
                                         RT = RT_reint[1], 
                                         mz = mz_reint[1]) , by = Feature ]
  f_rep[Feature %in% f_C57 , marker_C57 := T]
  f_rep[Feature %in% f_DBA , marker_DBA := T]
  f_rep[Feature %in% f_IS  , marker_IS  := T]
  f_rep[  , marker := (marker_C57 | marker_DBA | marker_IS) ,by = Feature]
  f_rep[Feature ==   rep_feature & marker == T, rep_marker := T ,by = Feature]
  
  # bind with metabolic cluster data
  dat <- merge(f_rep, clstr, by.x = "Feature", by.y = "id", all = T)
  
  # extrapolate metabolic correlation groups from representative features to all associated features  
  dat[ , metb_cor_grp_extrap := as.integer(mean(metb_cor_grp, na.rm = T)), by = rep_feature]
  dat[is.nan(metb_cor_grp_extrap) , metb_cor_grp_extrap := NA]
  col.del <- colnames(dat)[-1]
  
  # get id data
  dat.id <- data.table(read_xlsx(id.file, sheet = "ALL_FEATURES"))
  dat.id[ , c(col.del) := NULL] # delete columns that are already present in dat
  dat.id$new_feautures <- "no"  # set all features that are already present in previous id.file to "no"
  dat.id <- merge(dat.id, dat, by = "Feature", all.y = T,sort = F)
  dat.id[is.na(new_feautures), new_feautures := "yes"]
  
  # export new id data
  write.xlsx(dat.id,paste0(dat.loc, "Results/R_export/FEATURE_ID.xlsx"),
             showNA = F,
             row.names = F,
             sheetName = "ALL_FEATURES")
  
  #out
  return(dat.id)
  
  }