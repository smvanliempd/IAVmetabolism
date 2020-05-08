#  prepare data for stan models per metabolite
get.stan.data <- function(id.data, chem.cor.out) {
  
  # get relevant features and data table
  fts <- id.data[Use == T, Feature]
  dat <- chem.cor.out$data[Sample.Group != "QC" & is.na(Filler) & Feature %in% fts, 
                           .(Feature,
                             Sample.Group,
                             Subject.Species,
                             Challenge,
                             Animal,
                             Time.Factor,
                             Time.Day,
                             Weight,
                             Area_reint_mfc_qc) ]
  
  
  #  prepare data per metabolite
  dat.stan <- sapply(fts, function(f){
    
    d <- dat[Feature == f , .(Feature,Sample.Group,Animal,Subject.Species ,Challenge,Time.Factor,Time.Day,Weight,Area_reint_mfc_qc) ]
    d <- d[order(Subject.Species,-Challenge,Time.Factor)]
    
    # log trans
    d[ , Area_l:= log(Area_reint_mfc_qc)]
    
    # center/scale
    mu_data_l <- mean(d$Area_l)
    sd_data_l <- sd(d$Area_l)
    d[ , Area_lcs := (Area_l - mu_data_l)/sd_data_l]
    
    # set indices
    d[ , idx_animal  :=as.integer(factor(Animal, levels = unique(Animal)))]
    d[ , idx_strain  :=as.integer(factor(Subject.Species, levels = unique(Subject.Species))) ]
    
    # define all groups (n = 16)
    d[ , idx_group  := .GRP, by = c("Time.Factor","Subject.Species","Challenge")]
    
    # check groups
    N_gr <-  d[ , max(idx_group)]
    
    # prepare stan data
    dl <- list(N         = nrow(d),
               N_gr      = N_gr,
               mu_data_l = mu_data_l,
               sd_data_l = sd_data_l,
               idx_an    = d[ , idx_animal],
               idx_st    = d[ , idx_strain],
               idx_gr    = d[ , idx_group],
               time      = d[ , Time.Day],
               G         = sapply(1:N_gr, function(i) ifelse( d[ , idx_group] == i,1L,0L) ) ,
               Area      = d[ ,Area_reint_mfc_qc],
               Area_lcs  = d[ ,Area_lcs] )
    
    # out
    return( dl )
  },USE.NAMES = T,simplify = F)
  
  return(dat.stan)
  
}