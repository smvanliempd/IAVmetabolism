# cleaning, adjustment and modellinng of auto-integrated and manually integrated LCMS data
# ptm <- proc.time()
source("paths.R")
source("load_functions.R")

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
dat.merge <- merge.fun(dat.pos = dat.pos, dat.neg = dat.neg); rm(dat.pos,dat.neg) ; gc()
dat.merge <- fill.fun( qccorr.output = dat.merge )
dat.merge <- log.scale.fun(fill.out = dat.merge)

# model, post-hoc tests, feature selection
dat.mod <- mod.fun( stnd.out = dat.merge ); rm(dat.merge) ; gc()
dat.mod <- posthoc.fun( mod.out = dat.mod, pars.file = file.parm.all, contrast.file = file.contrast )
dat.mod <- select.fun( posthoc.out = dat.mod, pars.file = file.parm.all, contrast.file = file.contrast )

# merge reintegrated data
dat.reint <- reint.bind.fun(select.output = dat.mod, meta.file = file.meta.all, reint.files = file.reint)#; rm(dat.mod); gc()
dat.reint <- reint.clean.fun( reint.output = dat.reint, pars.file = file.parm.all)
dat.reint <- reint.adjust.fun( reint.clean.output = dat.reint, pars.file = file.parm.all)
dat.reint <- reint.log.scale.fun( reint.adjust.out = dat.reint)

# modeling reintegrated data
dat.mod.reint <- reint.mod.fun(stnd.out = dat.reint,pars.file = file.parm.all, contrast.file = file.contrast )
dat.mod.reint <- reint.select.fun( reint_ph.output = dat.mod.reint, pars.file = file.parm.all, contrast.file = file.contrast )

# correlation analysis
dat.corr <- chem.cor.fun(dat.merge.out = dat.mod.reint)
dat.corr <- metb.cor.fun(dat.corchem.out = dat.corr)

# Save the whole data thing to an .rds file
saveRDS(dat.corr, file = paste0(getwd(),"/IVAmetabolism_data.rds"))

# bind and export reintegrated feature properties and identifications
dat.id <- id.bind.fun(metb.cor.out = dat.corr, id.file = file.feature )

# 32.5 minutes on i7-9700K/3.60GHz/32Gb RAM
# proc.time() - ptm

# Bayesian multi-level models pre feature
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# prepare data for stan models
dat.stan <- get.stan.data(id.data = dat.id,chem.cor.out = dat.corr)

# compile model
mod.samples <- stan_model( file = paste0(getwd(),"/Models/msg_20190904_001_cc.stan" ))

# run models and extract samples
dat.samples <- get.stan.samples(dat = dat.stan, mod = mod.samples)

# Bind sample data to identifications
# this function can also be used for updating identifications
dat.samples.id <- id.bind.samples.fun(samples.out = dat.samples,meta.file = file.meta.all,id.file = file.feature )

# plotting
source("plotting.R")



