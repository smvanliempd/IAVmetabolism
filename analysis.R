
source("functions.R")
source("paths.R")

# Read and clean raw data for positive ionization mode
dat.pos <-  read.fun( file.path = file.raw.pos, meta.file = file.meta.all, polarity = "POS" )
dat.pos <-  clean.fun( read.output  = dat.pos,pars.file = file.parm.all, polarity = "POS" )
dat.pos <-  mfc.fun( clean.output = dat.pos )
dat.pos <-  qccorr.fun( mfc.output = dat.pos,pars.file = file.parm.all, polarity = "POS" )

# Read and clean raw data for negative ionization mode
dat.neg <-  read.fun( file.path  = file.raw.neg, meta.file = file.meta.all, polarity = "NEG" )
dat.neg <-  clean.fun( read.output = dat.neg, pars.file = file.parm.all, polarity = "NEG" )
dat.neg <-  mfc.fun( clean.output = dat.neg )
dat.neg <-  qccorr.fun( mfc.output  = dat.neg, pars.file = file.parm.all, polarity = "NEG" )

# merge polarities, fill and logtrans/standardize
dat.merge <- merge.fun(dat.pos = dat.pos, dat.neg = dat.neg); rm(dat.pos,dat.neg); gc()
dat.merge <- fill.fun( qccorr.output = dat.merge )
dat.merge <- log.scale.fun(fill.out = dat.merge)

# model, post-hoc tests, feature selection
dat.mod <- mod.fun( stnd.out = dat.merge, pars.file = file.parm.all )#; rm(dat.merge) ; gc()
dat.mod <- posthoc.fun( mod.out = dat.mod, contrast.file = file.contrast ) # ~25 min 
dat.mod <- select.fun( posthoc.out = dat.mod, pars.file = file.parm.all, contrast.file = file.contrast )
