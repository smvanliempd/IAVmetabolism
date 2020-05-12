# plot heatmaps for selcted metabolites 
plot.heatmap <- function(dat_plot, set, class_order ,cut_low, cut_high, plot_magnitude = T, dvc = "png") {
  
  # define variables for heatmaps
  vars.fc  <- c("del_pl_bl","del_d_3i_0i","del_c_3i_0i","del_c_5i_0i","del_c_8i_0i","del_c_18i_0i","del_c_30i_0i")
  vars.mag <- c("I_C_00") 
  
  # define fold-change and magnitude observations for heatmaps
  obs.fc  <- c("Final_ID","variable","strain","baseline","time",paste0("Class_",set),"P_pos","q50")
  obs.mag <- c("Final_ID","variable","group","strain","time",paste0("Class_",set),"P_pos","q50")
  
  # extract data
  dat.ftr <- dat_plot$features
  dat.fc  <- dat_plot$dat_deltas[ variable %in% vars.fc, .SD, .SDcols = obs.fc] # data for fold change
  dat.mag <- dat_plot$dat_groups[ group %in% vars.mag, .SD, .SDcols = obs.mag]  # data for absolute intensities for the t=0 C57BL/6J strain
  dat <- rbind(dat.fc, dat.mag, fill = T)
  
  # change name of Class column for faceting later
  setnames(dat,paste0("Class_",set),"Class")
  
  # delete rows for non-included mets
  dat <- dat[!is.na(Class)]
  
  # Normalize absolute intensities of the magnitude 
  dat[group == vars.mag, q50 := q50/median(q50), by = "Class" ]
  
  # get number of metabolites to plot
  n_mets <- length(unique(dat$Final_ID))
  
  # clean up and set levels
  dat[ , Class := factor(Class, class_order)]
  dat[ , Facet_vert := ifelse(variable == "del_d_3i_0i","DBA/2J",
                              ifelse(variable %in% c("del_c_3i_0i","del_c_5i_0i","del_c_8i_0i","del_c_18i_0i","del_c_30i_0i"),"C57BL/6J",
                                     ifelse(variable == "del_pl_bl","Inter-strain","Magnitude")))]
  dat[ , Facet_vert := factor(Facet_vert, c("Inter-strain","DBA/2J","C57BL/6J","Magnitude"))]
  dat[ , Time_discr := as.factor(time)]
  dat[ is.na(Time_discr), Time_discr := "BL"]
  dat[ Time_discr == "0", Time_discr := "M"]
  dat[ ,Time_discr := factor(Time_discr, levels = c("30","18","8","5","3","BL","M") )]
  
  # # convert log2 fold changes (L2FC) to percentages increase/decrease (percent change, P_delta) with respect to baseline (=100%)
  #   P_delta  = 2^L2FC - 1
  #   example: D2/B6 -> 50/200 = 0.25 -> L2FC =log2(0.25) = -2.0 so perc decrease is 100*(2^(0.25) - 1) = -75% from baseline
  #   example: D2/B6 -> 600/300 = 2.0 -> L2FC =log2(2.0)  =  1.0 so perc decrease is 100*(2^(1.0) - 1)  = +100% from baseline
  dat[ , q50_cens := 100*( ( 2^q50 ) - 1 )]
  dat[q50_cens <  cut_low, q50_cens  :=  cut_low]
  dat[q50_cens >  cut_high, q50_cens :=  cut_high]
  
  # define tile sizes, based on P_pos.
  # sqrt because area needs to map to confidence and not the height/width of the tile
  dat[P_pos > 0.5 , Tile_size :=  sqrt(2*(P_pos - 0.5)) ] #2*(P_pos - 0.5) ]#
  dat[P_pos < 0.5 , Tile_size :=  sqrt(1 - 2*P_pos) ] # 1 - 2*P_pos]#
  
  # set ID level according to P_pos
  # met.levels <- dat[variable == vars.fc[1] , .(P_pos = P_pos), by = Final_ID]
  # met.levels[ , met_levels := order(P_pos)]
  # met.levels <- met.levels[ , Final_ID[met_levels]]
  # dat[ ,Final_ID := factor(Final_ID, met.levels)]
  
  # set ID level to manual order
  met.levels <- dat.ftr[,Final_ID]
  dat[ ,Final_ID := factor(Final_ID, met.levels)]
  
  # plot
  if (plot_magnitude == F) {
    dat <- dat[Facet_vert != "Magnitude"]
  }
  p <- ggplot(dat,aes(y = Time_discr, x =  Final_ID, fill = q50_cens)) +
    geom_tile(aes(width = Tile_size, height = Tile_size), col = "grey80") +
    geom_hline(yintercept = (0:6)+0.5, size = 0.25, color = "grey90") +
    geom_vline(xintercept = (0:(n_mets))+0.5, size = 0.25, color = "grey90") +
    scale_fill_gradient2(name = expression(paste(Delta,"%")), 
                         low = "blue", mid = "white", high = "red",
                         limits = c(cut_low, cut_high),
                         guide=guide_colorbar(raster=F,nbin = 17) )+
    labs(y = "Days post-infection") +
    facet_grid(Facet_vert ~ Class , scales = "free", space = "free") +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          panel.grid = element_blank(),
          legend.key.height = unit(.15,"in"), 
          axis.text.x = element_text(angle = 90, vjust = .3, hjust = 1),
          legend.justification = c(0,1),
          strip.text.y = element_text(angle = 0),
          aspect.ratio = 1)
  ggsave(paste0(getwd(),"/Plots/Heatmaps/",set,"_heatmap_perc.", dvc), p,
         width = .15 * n_mets + 2,
         height = 5,
         units = "in", dpi = 600,device = dvc)
}
