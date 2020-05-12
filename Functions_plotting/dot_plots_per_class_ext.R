dot.plots.per.class.ext <-  function(dat_plot, set,  p_effect = 0.8, sel_effect = F, dvc = "png") {
  
  require(egg)
  require(grid)
  require(cowplot)
  
  # define variables for heatmaps
  vars.fc  <- c("del_pl_bl","del_d_3i_0i","del_c_3i_0i","del_del_str_inf") #,"del_d3i_c3i"
  
  # define observations for heatmaps
  obs.fc  <- c("Final_ID","variable","strain","baseline","time",paste0("Class_",set),"P_pos","q5","q25","q50","q75","q95")
  
  # extract data 
  dat     <- dat_plot$dat_deltas[ variable %in% vars.fc, .SD, .SDcols = obs.fc] # data for fold change
  dat.ftr <- dat_plot$features
  
  # change name of Class column for faceting later
  setnames(dat,paste0("Class_",set),"Class")
  
  # delete rows for non-included mets
  dat <- dat[!is.na(Class)]
  
  # get classes in set
  classes <- na.exclude(unique(dat$Class) )
  
  # plot per class
  pp <- sapply(classes , function(cls) {
    
    # get data per class
    dat.cls <- dat[Class == cls]
    
    # clean up and set levels
    dat.cls[ , Facet_vert := ifelse(variable == "del_pl_bl", "IS Healthy",ifelse(variable =="del_del_str_inf" , "IS Delta 3 dpi" ,"3 dpi"))]
    dat.cls[ , Facet_vert := factor(Facet_vert, levels = c("IS Healthy", "3 dpi","IS Delta 3 dpi"))]
    
    # set ID level according to P_pos
    # met.levels <- dat.cls[variable == vars.fc[1] , .(P_pos = P_pos), by = Final_ID]
    # met.levels[ , met_levels := order(P_pos)]
    # met.levels <- met.levels[ , Final_ID[met_levels]]
    # dat.cls[ ,Final_ID := factor(Final_ID, met.levels)]
    
    # set ID level to manual order
    met.levels <- dat.ftr[,Final_ID]
    dat.cls[ ,Final_ID := factor(Final_ID, met.levels)]
    
    # highlight effects with reliabale effect directions excluding "del_pl_bl"
    dat.cls[ , flag_direction := ifelse((P_pos <= p_effect | P_pos >= 1 - p_effect), "x", "o")] #variable %in% vars.fc[1] &
    fts <- droplevels(dat.cls[flag_direction == "x", Final_ID])
    
    # prepare plot data and parameters
    dw <- .5
    
    #function for changing log2(fc) in percentage delta change
    lfc2perc <- function(a) { 100*( ( 2^a ) - 1 ) } 
    if (sel_effect == T) dat.cls <- dat.cls[Final_ID %in% fts ]
    
    # get number of metabolites to plot
    n_mets <- length(unique(dat.cls$Final_ID))
    
    # plot
    p1 <- ggplot(dat.cls[variable != "del_del_str_inf"], aes(x= Final_ID, y = q50)) + 
      geom_hline(yintercept = 0, lty = 2, col = "red") +
      geom_linerange(aes(ymin = q25, ymax = q75, col = strain),
                     position = position_dodge(width = dw), size = 1 )+
      geom_linerange(aes(ymin = q5, ymax = q95, col = strain),
                     position = position_dodge(width = dw) , alpha = 0.7)+
      geom_point(aes( col = strain, shape = flag_direction),
                 position = position_dodge(width =dw), fill = "white")+
      scale_color_manual(name = NULL,values = c("blue","#ff8c00","black") ,
                         labels = c("C57 Infected","DBA Infected","IS Healthy"))+
      scale_shape_manual(values = c(21,19), guide = F)+
      labs( y = expression(paste(Delta,"% (50/90 PI)") ) ) +
      scale_y_continuous(breaks = c(-log2(c(2,5,10)),log2(c(1,2,5,11))),labels =  lfc2perc )+
      facet_grid(Facet_vert~Class,scales = "free", space = "free") +
      theme_bw() +
      theme(axis.text.x = element_blank(), axis.ticks.length.x = unit(0,"in"),
            axis.title.x = element_blank(), panel.grid.minor = element_blank())
    
    # set fixed panel sizes
    p1 <- set_panel_size(p1, width = unit(.2 * n_mets, "in"), height = unit(1.5, "in") )
    
    # save plot
    cls.name <-  str_replace(make.names(cls), "\\.", "_")
    ggsave2(paste0(getwd(),"/Plots/Dot_plots/Extended/dots_per_class_ext_",set,"_",cls.name,"_infected.", dvc), p1,
            width = .2 * n_mets + 5,
            height = 8,
            units = "in", dpi = 600,device = dvc)
    
    p2 <-  ggplot(dat.cls[variable == "del_del_str_inf"], aes(x= Final_ID, y = q50)) + 
      geom_hline(yintercept = 0, lty = 2, col = "red") +
      geom_linerange(aes(ymin = q25, ymax = q75, col = strain),
                     position = position_dodge(width = dw), size = 1 )+
      geom_linerange(aes(ymin = q5, ymax = q95, col = strain),
                     position = position_dodge(width = dw) , alpha = 0.7)+
      geom_point(aes( col = strain, shape = flag_direction),
                 position = position_dodge(width =dw), fill = "white")+
      scale_color_manual(name = NULL,values = c("red") ,
                         labels = c("IS Infected"))+
      scale_shape_manual(values = c(21,19), guide = F)+
      labs( y = expression(paste(Delta,Delta,"% (50/90 PI)") ) ) +
      facet_grid(Facet_vert~.,scales = "free", space = "free") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = .3, hjust = 1),
            axis.title.x = element_blank(), panel.grid.minor = element_blank() )
    
    p2 <- set_panel_size(p2, width = unit(.2 * n_mets, "in"), height = unit(1.5, "in") )
    
    ggsave2(paste0(getwd(),"/Plots/Dot_plots/Extended/dots_per_class_ext_",set,"_",cls.name,"_infdiff.", dvc), p2,
            width = .2 * n_mets + 5,
            height = 8,
            units = "in", dpi = 600,device = dvc)
    
    
  },simplify = F, USE.NAMES = T)
  
}
