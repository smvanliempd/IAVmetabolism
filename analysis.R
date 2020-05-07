# cleaning, adjustment and modellinng of auto-integrated and manually integrated LCMS data
source("paths.R")
sapply(list.files(fun.loc), function(f) source(paste0(fun.loc,f)))


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
dat.mod <- mod.fun( stnd.out = dat.merge, pars.file = file.parm.all ); rm(dat.merge) ; gc()
dat.mod <- posthoc.fun( mod.out = dat.mod, contrast.file = file.contrast )
dat.mod <- select.fun( posthoc.out = dat.mod, pars.file = file.parm.all, contrast.file = file.contrast )

# merge reintegrated data
dat.reint <- reint.bind.fun(select.output = dat.mod, meta.file = file.meta.all, reint.files = file.reint)#; rm(dat.mod); gc()
dat.reint <- reint.clean.fun( reint.output = dat.reint, pars.file = file.parm.all)
dat.reint <- reint.adjust.fun( reint.clean.output = dat.reint, pars.file = file.parm.all)
dat.reint <- reint.log.scale.fun( reint.adjust.out = dat.reint)

# modeling reintegrated data
dat.mod.reint <- reint.mod.fun(stnd.out = dat.reint,pars.file = file.parm.all, contrast.file = file.contrast )   # ~ 15 minutes
dat.mod.reint <- reint.select.fun( reint_ph.output = dat.mod.reint, pars.file = file.parm.all, contrast.file = file.contrast )

# correlation analysis
dat.corr <- chem.cor.fun(dat.merge.out = dat.mod.reint)
dat.corr <- metb.cor.fun(dat.corchem.out = dat.corr)

# Save the whole data thing to an .rds file
saveRDS(dat.corr, file = paste0(getwd(),"/IVAmetabolism_data.rds"))

# export reintegrated feaure properties and identifications
dat.export <- export.fun(dat.corr, "FEATURE_ID_v16_enriched_16.xlsx" )



