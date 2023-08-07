#' Create a custom publication theme for ggplot
#'
#' @param base_size A number to set the base font size. Default is 14.
#' @param base_family A character string specifying the base font family. Default is empty.
#' @return A complete ggplot theme with the specified characteristics.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Length, Petal.Length)) + geom_point() + theme_Publication()
#' }
#' @importFrom ggplot2 theme element_text element_rect element_line rel element_blank
#' @importFrom ggthemes theme_foundation
#' @importFrom grid unit

theme_Publication <- function(base_size=14, base_family="") {
  ggthemes::theme_foundation(base_size=base_size, base_family=base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.2, "cm"),
      #legend.margin = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin = grid::unit(c(10,5,5,5),"mm"),
      strip.background = element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    )
}
