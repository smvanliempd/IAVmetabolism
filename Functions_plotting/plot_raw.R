plot.raw <- function(dat, dvc = "png") {
  fts <- dat$features$Feature
  sapply(fts, function(f_raw) {
    d_plot <- dat$dat_groups[Feature == f_raw[1]]
    d_raw  <- dat$dat_raw[Feature == f_raw[1]]
    d_plot <- merge(d_plot, d_raw[,.(Sample.Group,Animal,Area_reint_mfc_qc)], by.x = "group", by.y = "Sample.Group" )
    m_name <- d_plot[Feature == f_raw, unique( Final_ID)]
    dw <- 1
    p <- ggplot(d_plot, aes( x = time, col = Challenge, fill = Challenge)) + #,group = Animal
      geom_point(aes( y = Area_reint_mfc_qc, pch = "r"), size = 1, alpha = 0.4, position = position_dodge(width = dw)) +
      geom_point(aes( y = q50, pch = "m"), position = position_dodge(width = dw)) +
      geom_line(aes( y = q50))+
      geom_linerange(aes(ymin = q5, ymax = q95), alpha = 0.05, size = 1, position = position_dodge(width = dw)) +
      geom_linerange(aes(ymin = q25, ymax = q75), alpha = 0.05, size = 2, position = position_dodge(width = dw)) +
      scale_color_manual(name = "Challenge", values = c(infection = "#ff0033", mock = "#000000")) + 
      scale_fill_manual(name = "Challenge", values = c(infection = "#ff0033", mock = "#000000")) +
      scale_x_continuous(breaks = c(0,3,5,8,18,30))+
      # scale_y_log10() +
      scale_shape_manual(values = c(17,1), name="Data type",label=c(expression(paste(mu["Q50"]," model")),"raw data") )+
      labs(title = paste0(m_name, " (", f_raw,")"), x = "Time Post-infection (Days)", 
           y = expression(paste("Signal [adjusted, PI50/90"[mu],"]"))) +
      facet_wrap(~strain, nrow = 2) +
      theme_bw() +
      theme(panel.grid.minor  = element_blank())
    
    m_name <- make.names(m_name)
    m_name <- str_replace_all(m_name,"\\.","_")
    ggsave2(filename = paste0(getwd(),"/Plots/Raw/",m_name,"_raw.",dvc),
            plot = p,width = 5,height = 4, units = "in", device = dvc, dpi = 600)
  } )
}
