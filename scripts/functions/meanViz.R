lWideToLong <- function(data, ...){
  reshape2::melt(
    data, 
    id.vars = c("PID", "TID", "TIDnum")
  )
}


statVarPlot <- function(dataIn, name, namList, stat = "mean", ...){
  # dataIn <- dtWorker$viz$MoodOverall
  # name <- "Mood"
  # namList <- varNameListWorker
  # stat <- "mean"
  
  varNames <- names(dataIn)[!names(dataIn) %in% c("PID", "TID", "TIDnum")]
  
  dataLong <- dataIn %>%
    mutate(TID = stri_replace_all_regex(
      TID,
      pattern = c('Morning', 'Afternoon'),
      replacement = c('12:00:00', '19:00:00'),
      vectorize = FALSE
    ) %>% as.POSIXct) %>%
    lWideToLong %>%
    mutate(variable = recode(variable, !!!namList))
  
  dataLongNaRm <- dataLong %>% filter(!is.na(value))
  
  nTimeLong <- dataLongNaRm %>%
    select(-PID) %>%
    group_by(TIDnum, TID, variable) %>%
    summarise(
      n = n()
    ) %>%
    mutate(
      perc = n/length(unique(dataLongNaRm$PID))*100
    ) 
  
  
  
  pTidTop <- ggplot(nTimeLong, aes(x = TIDnum, y = perc, color = variable, linetype = variable)) +
    geom_line(
      stat="smooth",
      method = "loess",
      formula = 'y ~ x',
      span = 0.1,
      se = FALSE,
      #size = 1,
      #color = "black",
      alpha = 0.75
    ) +
    ggthemes::scale_colour_calc() +
    geom_hline(yintercept = 80,
               linetype = "dashed",
               color = "black") +
    annotate(
      # add white background
      "label",
      x = max(dataLong$TIDnum),
      y = 80,
      label = "80%",
      vjust = 0.5,
      hjust = -1,
      label.size = NA
    )  +
    coord_cartesian(clip = "off") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(y = "Percentage") +
    #scale_colour_manual(values = rep("black", 20)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 10, 0, 10, "mm")
    )
  pTidMain <-
    ggplot(dataLongNaRm,
           aes(x = TIDnum, y = value, color = variable)) +
    #geom_point(alpha = 0) +
    stat_summary(fun = stat, geom = "line") +
    labs(y = paste0(stat),
         x = "Time Index") +
    #scale_colour_manual(values = RColorBrewer::brewer.pal(length(unique(dataLong$variable)), "Set3")) +
    ggthemes::scale_colour_calc() +
    scale_x_continuous(breaks = seq(0, max(dataLongNaRm$TIDnum), 25)) +
    scale_y_continuous(limits = switch((stat == "mean")+1, NULL, c(min(dataLongNaRm$value), max(dataLongNaRm$value))),
                       breaks = switch((stat == "mean")+1, waiver(), seq(min(dataLongNaRm$value), max(dataLongNaRm$value), (max(dataLongNaRm$value)-min(dataLongNaRm$value))/5))) +
    theme_Publication() +
    theme(
      panel.border = element_rect(colour = "black"),
      plot.margin = margin(0, 10, 10, 10, "mm")
    )
  pTidTitle <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste("Variable Set:", name, "[Time Id Plot]"),
      fontface = 'bold',
      x = 0.5,
      hjust = 0.5
    )
  
  pTid <- cowplot::plot_grid(
    pTidTitle,
    pTidTop,
    pTidMain,
    align = "v",
    nrow = 3,
    rel_heights = c(0.1, 1 / 5, 4 / 5)
  )
  
  
  
  pDateTop <- ggplot(nTimeLong, aes(x = TID, y = perc, color = variable, linetype = variable)) +
    geom_line(
      stat="smooth",
      method = "loess",
      formula = 'y ~ x',
      span = 0.1,
      se = FALSE,
      #size = 1,
      #color = "black",
      alpha = 0.75
    ) +
    ggthemes::scale_colour_calc() +
    geom_hline(yintercept = 80,
               linetype = "dashed",
               color = "black") +
    annotate(
      # add white background
      "label",
      x = max(dataLong$TID),
      y = 80,
      label = "80%",
      vjust = 0.5,
      hjust = -1,
      label.size = NA
    )  +
    coord_cartesian(clip = "off") +
    scale_y_continuous(limits = c(0, 100)) +
    labs(y = "Percentage") +
    #scale_colour_manual(values = rep("black", 20)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 10, 0, 10, "mm")
    )
  pDateMain <-
    ggplot(dataLongNaRm,
           aes(x = TID, y = value, color = variable)) +
    #geom_point(alpha = 0) +
    stat_summary(fun = stat, geom = "line") +
    labs(y = paste0(stat),
         x = "Date") +
    ggthemes::scale_colour_calc() +
    scale_x_datetime(breaks = scales::date_breaks("10 day")) +
    scale_y_continuous(limits = switch((stat == "mean")+1, NULL, c(min(dataLongNaRm$value), max(dataLongNaRm$value))),
                       breaks = switch((stat == "mean")+1, waiver(), seq(min(dataLongNaRm$value), max(dataLongNaRm$value), (max(dataLongNaRm$value)-min(dataLongNaRm$value))/5))) +
    theme_Publication() +
    theme(
      panel.border = element_rect(colour = "black"),
      plot.margin = margin(0, 10, 10, 10, "mm")
    )
  pDateTitle <- cowplot::ggdraw() +
    cowplot::draw_label(
      paste("Variable Set: Key-Need", "[Date Plot]"),
      fontface = 'bold',
      x = 0.5,
      hjust = 0.5
    )
  
  pDate <- cowplot::plot_grid(
    pDateTitle,
    pDateTop,
    pDateMain,
    align = "v",
    nrow = 3,
    rel_heights = c(0.1, 1 / 5, 4 / 5)
  )
  
  list(pTid, pDate)
  pTid
}
