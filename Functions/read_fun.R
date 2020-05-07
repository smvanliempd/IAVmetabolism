# Read the data from a MarkerLynx/EZInfo output .txt file
read.fun   <- function( file.path, meta.file, polarity) {
  
  #read raw data
  dat.raw  <- read.csv(file.path, sep = "\t", stringsAsFactors = F)
  
  # read meta data
  dat.meta <-  read_xlsx( meta.file)
  colnames(dat.meta) <- make.names(colnames(dat.meta))
  dat.meta <- subset(dat.meta, Polarity == polarity & is.na(Exclude) )
  
  # select all columns containing features
  dat.feat <- grep("[X,N,P][0-9]{1}", colnames(dat.raw))
  
  # Separate rtmz values and intensity values,
  # delete previous meta data and rename 1st column
  dat.rtmz  <- dat.raw[  1:2  , dat.feat]
  dat.clean <- dat.raw[-(1:2) , c(1,dat.feat)]
  colnames(dat.clean)[1] <- "File.Name"
  
  # bind values with new meta data
  dat.full <- left_join(dat.meta , dat.clean, by = "File.Name")
  
  # select all columns containing features
  dat.feat <- grep("[X,N,P][0-9]{1}", colnames(dat.full))
  
  # extract feature names, replace with accurate polarity identifier and update column names
  feat.names <- colnames(dat.full[ , dat.feat])
  feat.names <- str_replace(feat.names , "X", ifelse(polarity == "POS", "P", "N") )
  colnames(dat.full)[dat.feat] <- feat.names
  colnames(dat.rtmz) <- feat.names
  
  # transform to long format in data.table 
  dat.long <- data.table(gather(dat.full,key = "Feature", value = "Signal", dat.feat ))
  setkey(dat.long, Feature)
  
  # bind long data with rtmz data
  dat.rtmz <- data.table(Feature = colnames(dat.rtmz),
                         RT      = unlist(c(dat.rtmz[1,])),
                         mz      = unlist(c(dat.rtmz[2,])) )
  setkey(dat.rtmz, Feature)
  dat.long <- dat.long[dat.rtmz]
  
  return(dat.long)
}