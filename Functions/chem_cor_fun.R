# Check "chemical" correlations i.e. check for adducts and fragments
chem.cor.fun <- function( dat.merge.out ){
  
  # Function to check for adducts and fragments etc.
  # Not perfect but close enough:
  #  e.g. problems with features c("P0.66_226.1810","P0.66_296.1325", "P0.66_250.1991")
  #  due to lockmass ineterference and low inetensities
  
  # cross check function
  cross_check <- function( dat ) {
    
    # this function checks for features that are correlated via another feature
    
    # initalization
    l  <- max(dat$chem_cor_grp_1)
    g2 <- sapply(1:l , function(i)  dat[chem_cor_grp_1 == i, chem_cor_grp_2])
    g1 <- sapply(g2  , function(g2) dat[chem_cor_grp_2 %in% g2, unique(chem_cor_grp_1)])
    cc <- list(g2 = g2, g1 = g1 )
    
    # loop cross-checks until there are no more changes in correlation groups
    # output (cc) are lists of correlated groups (chemically correlated)
    cc.recursive <- function(dat, cc ) {
      cc$g2p <- sapply(cc$g1  , function(g1) dat[chem_cor_grp_1 %in% g1, unique(chem_cor_grp_2)])
      cc$g1p <- sapply(cc$g2  , function(g2) dat[chem_cor_grp_2 %in% g2, unique(chem_cor_grp_1)])
      cc$g2  <- sapply(cc$g1p , function(g1) dat[chem_cor_grp_1 %in% g1, unique(chem_cor_grp_2)])
      cc$g1  <- sapply(cc$g2p , function(g2) dat[chem_cor_grp_2 %in% g2, unique(chem_cor_grp_1)])
      if (identical(cc$g1, cc$g1p) & identical(cc$g2, cc$g2p)) {
        return(cc) 
      } else {
        cc.recursive(dat, cc) 
      }
    }
    
    # list correlated features
    cc <- cc.recursive(dat, cc)
    f_cor <- sapply(1:length(g1), function(i) dat[chem_cor_grp_1 %in% cc$g1[[i]] | chem_cor_grp_2 %in% cc$g2[[i]],
                                                  sort(unique(c(Feature_1,Feature_2))) ] )
    return(f_cor)
    
  }
  
  # prepare data 
  f_reint <- dat.merge.out$feature_data$reint_incl
  dat.cor <- dat.merge.out$data
  
  # Select reintegrated features and transform if logged
  d_cor <- dat.cor[Feature %in%  f_reint & Sample.Group != "QC", .(Feature,RT_reint,mz_reint,Sample.ID, Sample.Group,Area_reint_mfc_qc) ]
  setkey(d_cor, RT_reint)
  
  # make correlation matrix
  d_cor <- dcast(d_cor, Sample.ID ~ Feature,value.var = "Area_reint_mfc_qc")
  d_cor <- as.matrix(d_cor[,- 1])
  d_cor <- cor(d_cor, use = "complete.obs")
  d_cor <- data.table(Feature = colnames(d_cor), d_cor)
  d_cor <- melt(d_cor, measure.vars  = d_cor$Feature)
  setnames(d_cor, c("variable","value"), c("Feature_2","correlation"))
  
  # long format and bind with rtmz data
  d_rt <- dat.cor[ ,unique(RT_reint), by = Feature]
  d_cor <- merge(d_cor, d_rt, by = "Feature", variable.factor = F)
  setnames(d_cor, c("V1"), c("RT_reint") )
  d_cor <- merge(d_cor, d_rt, by.x = "Feature_2", by.y = "Feature", variable.factor = F)
  setnames(d_cor, c("Feature_2","Feature","V1"), c("Feature_1","Feature_2","RT_reint_2") )
  
  # determine RT differences between corelation pairs and select all chemically related features
  dRT     <- 0.005 # maximum RT difference between correlation pairs for them to be regarded as chemically related
  min_cor <- 0.8
  
  d_cor[ , delta_RT := abs(RT_reint - RT_reint_2)]
  d_cor <- d_cor[ delta_RT < dRT & correlation > min_cor  ]  # you can tune these parameters
  
  # make groups by features --> this can go inside the cross_check() function
  setkey(d_cor, "Feature_1")
  d_cor[ , chem_cor_grp_1 := .GRP, by = Feature_1]
  setkey(d_cor, "Feature_2")
  d_cor[ , chem_cor_grp_2 := .GRP, by = Feature_2]
  
  # determine all correlated features
  f_cor <- cross_check( d_cor )
  d_cor[ , chem_cor_grp := NA]
  for (i in 1: length(f_cor) ) d_cor[ , chem_cor_grp := ifelse(Feature_1 %in% f_cor[[i]] | Feature_2 %in% f_cor[[i]]  , i, chem_cor_grp) ]
  d_cor[ !is.na(chem_cor_grp), chem_cor_grp := .GRP, by = chem_cor_grp]
  
  # annotate chem_cor features in total data
  dat.cor[ , chem_cor_grp := NA]
  for (i in 1: length(f_cor) ) dat.cor[ , chem_cor_grp := ifelse(Feature %in% f_cor[[i]], i, chem_cor_grp) ]
  dat.cor[ !is.na(chem_cor_grp), chem_cor_grp := .GRP, by = chem_cor_grp]
  
  # Find represenative features for each chem_cor group (feature with highest median intensity)
  dat.cor[ , med_int := median(reint_log_stand, na.rm = T), by = c("chem_cor_grp", "Feature")]
  dat.cor[ , rep_feature := unique(Feature[which.max(med_int)])[1] , by = chem_cor_grp]
  
  # get representative features
  f_rep  <- dat.cor[!is.na(rep_feature), unique(rep_feature) ]
  
  # get representative markers and bind with output
  f_C57 <- dat.merge.out$selections_reint$feat$select_C57
  f_DBA <- dat.merge.out$selections_reint$feat$select_DBA
  f_int <- dat.merge.out$selections_reint$feat$select_InterSpecies
  f_tot <- dat.merge.out$selections_reint$feat$select_total
  
  dat.merge.out$selections_reint$feat$select_C57_rep <- f_C57[f_C57 %in% f_rep]
  dat.merge.out$selections_reint$feat$select_DBA_rep <- f_DBA[f_DBA %in% f_rep]
  dat.merge.out$selections_reint$feat$select_InterSpecies_rep <- f_int[f_int %in% f_rep]
  dat.merge.out$selections_reint$feat$select_total_rep <- f_tot[f_tot %in% f_rep]
  
  # annotate chem_cor features in correlation matrix
  d_cor[ ,chem_cor_grp := NA ]
  for (i in 1: length(f_cor) ) d_cor[ ,chem_cor_grp := ifelse(Feature_1 %in% f_cor[[i]], i, chem_cor_grp) ]
  d_cor[ , chem_cor_grp := .GRP, by = chem_cor_grp]
  
  # return prep
  dat.merge.out$data <- dat.cor
  dat.merge.out$correlations$chemical <- d_cor
  dat.merge.out$feature_data$included$representative_features <- f_rep
  
  return(dat.merge.out)
  
} 