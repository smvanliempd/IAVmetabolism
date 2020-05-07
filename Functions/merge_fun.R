# merge positive and negative data
merge.fun <- function( dat.pos, dat.neg ) {
  
  dpos <- dat.pos$data
  dneg <- dat.neg$data
  
  # add ion mode
  dpos$mode <- "POS"
  dneg$mode <- "NEG"
  
  # mereg pos and neg
  dmerge <- rbind(dpos,dneg)
  
  # out
  out <- list(data = dmerge,
              features = unique(dmerge$Feature),
              meta = list(pos = dat.pos$feature_data, 
                          neg = dat.neg$feature_data))
  
  return (out)
  
}