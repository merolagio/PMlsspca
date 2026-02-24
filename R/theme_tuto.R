
#' @importFrom ggplot2 ggplot aes theme theme_minimal 
#' @importFrom ggplot2 theme_light theme_bw theme_classic
#' @importFrom ggplot2 geom_bar geom_tile geom_rect geom_line geom_point
#' @importFrom ggplot2 geom_abline geom_vline geom_text geom_smooth
#' @importFrom ggplot2 element_rect  element_blank element_text unit labs
#' @importFrom ggplot2 scale_y_discrete scale_y_continuous scale_x_discrete scale_x_continuous 
#' @importFrom ggplot2 scale_fill_manual scale_fill_grey scale_fill_gradientn xlab ylab
#' @importFrom ggplot2 facet_wrap facet_grid vars
#' @importFrom ggplot2 annotate coord_polar coord_flip 
#' @noRd
NULL


## ---- default ggplot theme ----
theme_tuto <- function(base_size = 12, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title.position = "plot",
      legend.position = "bottom",
      panel.background = element_rect(colour = "black")
    )
}


##scale_y_continuous(labels = scales::percent)
# make_colours(7, "p")
# make_colours(7, "cbb")
# make_colours(7, "gg")

# pl <- pl + ggplot2::scale_fill_manual(values = make_colours(n))

make_colours = function(n, pal = c("cbb", "ggplot")){
  if (grepl("^c", pal[1])){
    cbb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    if (n <= 8) 
      cols = cbb[1:n]
    else{
      warnings("cannot use bcc pallette with more than 8 colours. Switch to ggplot")
      pal = "ggplot"
    }
  }
  if(grepl("^p", pal[1])){
    if (n < 10){
      pf.pal <- rev(c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#EF6548", "#D7301F", "#B30000", "#7F0000"))
      cols = pf.pal[1:n]
    }
    else{
      warnings("cannot use pf pallette with more than 8 colours. Switch to ggplot")
      pal = "ggplot"
    }
  }
  if (grepl("^g", pal[1])){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols <- gg_color_hue(n)
  }
  return(cols)
}

