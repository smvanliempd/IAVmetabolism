# Check if two features are correlated (also via intermediate features)
metb.cor.fun <- function( dat.corchem.out ) {
  # metabolic correlations
  
  # cross check function
  cross_check <- function( dat ) {
    
    # this function checks for features that are correlated via another feature
    
    # assign numbers to metabolic clusters
    setkey(dat, "Feature_1")
    dat[ , metb_cor_grp_1 := .GRP, by = Feature_1]
    
    setkey(dat, "Feature_2")
    dat[ , metb_cor_grp_2 := .GRP, by = Feature_2]
    
    l  <- max(dat$metb_cor_grp_1)
    g2  <- sapply(1:l , function(i)  dat[metb_cor_grp_1 == i, metb_cor_grp_2] ) # get corresponding group numbers from group 2
    g1  <- sapply(g2  , function(g2) dat[metb_cor_grp_2 %in% g2, unique(metb_cor_grp_1)]) 
    cc <- list(g2 = g2, g1 = g1 )
    
    # loop cross-checks until there are no more changes in correlation groups
    # output (cc) are lists of correlated groups (chemically correlated)
    cc.recursive <- function(dat, cc ) {
      cc$g2p <- sapply(cc$g1  , function(g1) dat[metb_cor_grp_1 %in% g1, unique(metb_cor_grp_2)])
      cc$g1p <- sapply(cc$g2  , function(g2) dat[metb_cor_grp_2 %in% g2, unique(metb_cor_grp_1)])
      cc$g2  <- sapply(cc$g1p , function(g1) dat[metb_cor_grp_1 %in% g1, unique(metb_cor_grp_2)])
      cc$g1  <- sapply(cc$g2p , function(g2) dat[metb_cor_grp_2 %in% g2, unique(metb_cor_grp_1)])
      if (identical(cc$g1, cc$g1p) & identical(cc$g2, cc$g2p)) {
        return(cc) 
      } else { 
        cc.recursive(dat, cc) 
      }
    } 
    
    # list correlated features
    cc <- cc.recursive(dat, cc)
    f_cor <- sapply(1:length(g1), function(i) dat[metb_cor_grp_1 %in% cc$g1[[i]] | metb_cor_grp_2 %in% cc$g2[[i]],
                                                  sort(unique(c(Feature_1,Feature_2))) ] )
    return(f_cor)
    
  } 
  
  # prepare data 
  state <- c("baseline","infected")
  species <- c("C57BL/6") #, "DBA/2" omitted because not enough data points to get meaningfull correlations
  dat.cor <- dat.corchem.out$data
  dat.cor <- dat.cor[Feature == rep_feature ] #get only representative features of reintegrated data
  dat.cor[ , State := ifelse( Challenge == "mock" | Time.Day == 0, state[1] , state[2])  ]
  
  
  dat.net <- sapply( species, function(sp){
    sapply( state, function(st) {
      
      # Select reintegrated features and log transform
      dat1 <- dat.cor[ State == st & Subject.Species == sp, .(Feature,Sample.ID,Area_reint_mfc_qc) ] #, log_trans_reint
      dat1[ , Area_reint_mfc_qc := log10(Area_reint_mfc_qc)] #log_trans_reint == F
      
      # make correlation matrix
      dat1 <- dcast(dat1, Sample.ID ~ Feature,value.var = "Area_reint_mfc_qc")
      dat1 <- as.matrix(dat1[,- 1])
      dat1 <- cor(dat1, use = "complete.obs", method = "pearson") #"spearman"
      dat1[lower.tri(dat1)] <- NA
      dat1 <- data.table(Feature = colnames(dat1), dat1)
      dat1 <- melt(dat1, measure.vars  = dat1$Feature, variable.factor = F)
      setnames(dat1, c("Feature","variable","value"), c("Feature_1","Feature_2","correlation_A"))
      dat1 <- dat1[!is.na(correlation_A)]
      
      # return(dat1)
      
      # determine RT differences between corelation pairs 
      cor_neg <- -0.90
      cor_pos <-  0.90
      dat1 <- dat1[  (correlation_A < cor_neg | correlation_A > cor_pos) & correlation_A != 1 ]  #delta_RT > 0.06 & 
      
      # determine all correlated features
      f_cor <- cross_check( dat1 )
      dat1[ , metb_cor_grp := NA]
      for (i in 1: length(f_cor) ) dat1[ , metb_cor_grp := ifelse(Feature_1 %in% f_cor[[i]] | Feature_2 %in% f_cor[[i]]  , i, metb_cor_grp) ]
      dat1[ !is.na(metb_cor_grp), metb_cor_grp := .GRP, by = metb_cor_grp]
      
      # check the same features in other state
      dat2 <- dat.cor[State == ifelse(st == state[1], state[2], state[1]) & Subject.Species == sp, .(Feature,Sample.ID,Area_reint_mfc_qc) ] #, log_trans_reint
      dat2[ , Area_reint_mfc_qc := log10(Area_reint_mfc_qc)] #log_trans_reint == F
      
      # make correlation matrix
      dat2 <- dcast(dat2, Sample.ID ~ Feature,value.var = "Area_reint_mfc_qc")
      dat2 <- as.matrix(dat2[,- 1])
      dat2 <- cor(dat2, use = "complete.obs", method = "pearson" ) #"spearman"
      dat2[lower.tri(dat2)] <- NA 
      dat2 <- data.table(Feature = colnames(dat2), dat2)
      dat2 <- melt(dat2, measure.vars  = dat2$Feature, variable.factor = F)
      setnames(dat2, c("Feature","variable","value"), c("Feature_1","Feature_2","correlation_B"))
      dat2 <- dat2[!is.na(correlation_B)]
      
      #merge 
      dat2 <- merge(dat1,dat2,  by = c("Feature_1","Feature_2"), all.x = T) #[ ,.(Feature_1,Feature_2,correlation_B)]
      dat2 <- dat2[ !(correlation_B < cor_neg | correlation_B > cor_pos), correlation_B := NA]  #, correlation_B := NA
      
      # make edges for network plot
      setnames(dat2, c("Feature_1","Feature_2"), c("from","to"))
      dat2[ , color := ifelse(correlation_A < 0, "red", "blue") ]
      dat2[ , value := abs(correlation_A)]
      
      # make nodes for network plot
      f_cor    <- dat2[ , unique(c(from,to) ) ]
      f_remain <- dat2[!is.na(correlation_B) , unique(c(from,to) ) ] #also correlated in other state
      nodes <- data.table(id = f_cor)
      nodes[, title := id]
      nodes[id %in% f_remain, color.background := "lightgreen"]
      
      # make node labels based on metabolic groups
      labs1 <- dat2[ , .(from, metb_cor_grp)]
      labs2 <- dat2[ , .(to, metb_cor_grp)]
      
      # return(list(labs1 = labs1, labs2 = labs2, dat2 = dat2))
      
      labs <- merge(labs1, labs2, by.x = c("from","metb_cor_grp"), by.y = c("to","metb_cor_grp"), all = T , allow.cartesian = T)
      
      # return(labs)
      
      labs[ , n := 1:.N, by = from]
      labs <- labs[ n == 1]
      labs[ ,n := NULL]
      
      nodes <- merge(nodes, labs, by.x = "id", by.y = "from")
      nodes[ , label := paste0(id," (clstr_", metb_cor_grp,")")]
      
      list(edges = dat2, nodes = nodes)
      
    }, USE.NAMES = T, simplify = F )
  }, USE.NAMES = T, simplify = F  )
  
  # out
  dat.corchem.out$correlations$metabolic <- dat.net
  return(dat.corchem.out)
}