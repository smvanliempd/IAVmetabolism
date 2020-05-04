source("paths.R")
source("functions.R")

#read and clean raw data for positive ionization
dat.pos <-  read.fun( file.path = file.raw.pos, meta.file = file.meta.all, polarity = "POS" )
dat.pos <-  clean.fun( read.output  = dat.pos,pars.file = file.parm.all, polarity = "POS" )
dat.pos <-  mfc.fun( clean.output = dat.pos )
dat.pos <-  qccorr.fun( mfc.output = dat.pos,pars.file = file.parm.all, polarity = "POS" )


