plot.diff.strain.3dpi <- function(dat_plot, set, class_order, dvc = "png", scale  = 2) {
  
  # load ggforce package
  require(ggforce)
  
  # define variables for heatmaps
  vars.fc  <- c("del_d_3i_0i","del_c_3i_0i")
  
  # define observations for heatmaps
  obs.fc  <- c("Final_ID","variable","strain","baseline","time",paste0("Class_",set),"q25","q50","q75")
  
  # extract data 
  dat  <- dat_plot$dat_deltas[ variable %in% vars.fc, .SD, .SDcols = obs.fc] # data for fold change
  
  # change name of Class column for faceting later
  setnames(dat,paste0("Class_",set),"Class")
  
  # delete rows for non-included mets
  dat <- dat[!is.na(Class)]
  
  # spread quantile data over strains
  dat[,variable := NULL]
  dat <- dcast(dat, ...~ strain, value.var = c("q25","q50","q75"))
  
  # clean up and set levels
  dat[ , Class := factor(Class, class_order)]
  
  # Set diameter of ellipse based on the sign of the effect since ellipse only allows symetric axes.
  #   This is not really a problem since intervals are very symetric anyway.
  dat[ , diam_C57 := ifelse(q50_C57 < 0, -1*(q25_C57 - q50_C57), (q75_C57 - q50_C57)) ]
  dat[ , diam_DBA := ifelse(q50_DBA < 0, -1*(q25_DBA - q50_DBA), (q75_DBA - q50_DBA)) ]
  
  # set colors for geom_point based in the quadrants that they are in
  dat[ , p_col := ifelse(sign(q50_C57) == sign(q50_DBA), "black", "red") ]
  
  # plot ellipses and dots
  p <- ggplot(dat, aes(x = q50_C57, y = q50_DBA) ) +
    geom_ellipse(aes(x0 = q50_C57, y0 = q50_DBA, a = diam_C57, b = diam_DBA, angle = 0 ),fill = "black", alpha = 0.05, col = NA) +
    geom_hline(yintercept = 0, size = 0.25) +
    geom_vline(xintercept = 0, size = 0.25) +
    geom_hline(yintercept = c(-1,1),lty = 2, col = "grey80") +
    geom_vline(xintercept = c(-1,1), lty = 2,col = "grey80") +
    geom_abline(intercept = 0,slope = 1,col = "grey80")+
    geom_point(aes(col = p_col),size = 2, show.legend = F ) +
    geom_smooth(method = "lm", col= "red", size = 0.5, se = F)+
    scale_color_manual(values = c("black","red")) +
    scale_y_continuous(breaks = -10:10)+
    scale_x_continuous(breaks = -10:10)+
    labs( x = expression(paste(log[2], "(FC) in C57 at 3 dpi")),
          y =expression(paste( log[2], "(FC) in DBA at 3 dpi") ))+
    facet_wrap(Class ~ ., nrow = 1 ) +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text = element_text(size = 14))
  
  
  # calculate plot dimensions
  n_class <- length(class_order)
  delx <- diff(range(dat$q50_C57))
  dely <- diff(range(dat$q50_DBA))
  
  s <- scale
  w <-  s * n_class
  h <- (s * dely/delx ) + 0.1  * n_class
  
  # save
  ggsave(paste0(getwd(),"/Plots/Strain_diff/",set,"_diff_strain_3dpi.", dvc), p,
         width = w,
         height = h,
         units = "in", dpi = 600,device = dvc)
}
