
sdVarPlot <- function(dataIn, name, namList,...){
  # dataIn <- dtMedical$viz$KeyNeedFulfillment
  # name <- "Core Need Fulfillment"
  # namList <- varNameListMedical
  
  varNames <- names(dataIn)[!names(dataIn) %in% c("PID", "TID", "TIDnum")]
  
  dataIn %>%
    mutate(TID = stri_replace_all_regex(
      TID,
      pattern = c('Morning', 'Afternoon'),
      replacement = c('12:00:00', '19:00:00'),
      vectorize = FALSE
    ) %>% as.POSIXct) %>%
    select(-PID) %>%
    #filter_at(vars(varNames), all_vars(!is.na(.))) %>%
    group_by(
      TIDnum, TID
    ) %>%
    summarise(
      n = n(),
      #across(varNames, ~mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
      across(varNames, ~n(.x), .names = "n_{.col}")
    ) 

  
  
  
  pTidMain <-
    ggplot(dataLong,
           aes(x = TIDnum, y = value, color = variable)) +
    geom_point(alpha = 0) +
    stat_summary(fun = mean, geom = "line") +
    labs(y = "Mean per Day",
         x = "Time Index") +
    #scale_colour_manual(values = RColorBrewer::brewer.pal(length(unique(dataLong$variable)), "Set3")) +
    ggthemes::scale_colour_calc() +
    scale_x_continuous(breaks = seq(0, max(dataLong$TIDnum), 25)) +
    theme_Publication() +
    theme(
      panel.border = element_rect(colour = "black"),
      plot.margin = margin(0, 10, 10, 10, "mm")
    )
  pTidTop <-
    ggplot(dataLong %>% filter(!is.na(value)), aes(x = TIDnum, y = (..count.. / max(count)) * 100),
           color = variable) +
    geom_freqpoly(binwidth = 2) +
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
  
  pDateMain <-
    ggplot(dataLong,
           aes(x = TID, y = value, color = variable)) +
    geom_point(alpha = 0) +
    stat_summary(fun = mean, geom = "line") +
    labs(y = "Mean per Day",
         x = "Date") +
    #scale_colour_manual(values = RColorBrewer::brewer.pal(length(unique(dataLong$variable)), "Set3")) +
    ggthemes::scale_colour_calc() +
    scale_x_datetime(breaks = scales::date_breaks("10 day")) +
    theme_Publication() +
    theme(
      panel.border = element_rect(colour = "black"),
      plot.margin = margin(0, 10, 10, 10, "mm")
    )
  pDateTop <-
    ggplot(dataLong,
           aes(
             x = TID,
             y = (..count.. / max(count)) * 100,
             color = variable
           )) +
    geom_freqpoly(binwidth = 1 * 3600 * 24) +
    geom_hline(yintercept = 80,
               linetype = "dashed",
               color = "black") +
    annotate(
      # add white background
      "label",
      x = max(unique(dataLong$TID)),
      y = 80,
      label = "80%",
      vjust = 0.5,
      hjust = -1,
      label.size = NA
    )  +
    coord_cartesian(clip = "off") +
    labs(y = "Percentage") +
    scale_colour_manual(values = rep("black", 20)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(0, 10, 0, 10, "mm")
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
