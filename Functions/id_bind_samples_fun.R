# bind samples with identifications and split in deltas and sample groups
id.bind.samples.fun <- function(samples.out, id.file , meta.file) {
  
  # get meta data for deltas and sample groups
  dm_del <- data.table(read_xlsx( meta.file, sheet = "DELTA"))
  dm_grp <- data.table(read_xlsx( meta.file, sheet = "GROUP"))
  dm_fts <- data.table(read_xlsx( id.file, sheet = "ALL_FEATURES"))
  
  # get sample quantile data
  dat <- samples.out$s_quantiles
  
  # bind quantile data with id data
  dat <- merge(dat, dm_fts, by = "Feature", all.x = T)
  
  # split quantile data in deltas and sample group
  d_del <- merge(dat, dm_del,  by.x = "variable", by.y = "delta.pars")
  d_grp <- merge(dat, dm_grp,  by.x = "variable", by.y = "group.pars")
  
  dat.out <- list(deltas = d_del, groups = d_grp)
  
  return(dat.out)
  
}