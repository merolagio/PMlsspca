## .Rprofile  (place in the project root)

local({
  ## Quiet startup
  # options(
  #   repos = c(CRAN = "https://cloud.r-project.org"),
  #   stringsAsFactors = FALSE
  # )
  
  ## ---- Packages ----
  suppressPackageStartupMessages({
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      library(ggplot2)
      
      ## ---- Your default ggplot theme ----
      ## Edit this function to match your style.
      theme_tuto <- function(base_size = 12, base_family = "") {
        ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            plot.title.position = "plot",
            legend.position = "bottom",
            panel.background = element_rect(colour = "black")
          )
      }

      
    } else {
      message("Package 'ggplot2' not installed. Run: install.packages('ggplot2')")
    }
  })

})
