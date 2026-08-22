suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(readxl)
  library(ggpubr)
  library(ggsci)
  library(gridExtra)
  library(grid)
  library(scales)
})

set.seed(3)

options(
  vsc.dev.args = list(width = 1500, height = 1500, pointsize = 12, res = 300),
  ggplot2.discrete.colour = c("#00AFBB", "#E7B800", "#FC4E07", "grey"),
  ggplot2.discrete.fill = c("#00AFBB", "#E7B800", "#FC4E07", "grey")
)

theme_set(theme_bw(base_size = 14) + theme(legend.title = element_blank()))

`%notin%` <- Negate(`%in%`)

drug_list <- c("CFZ", "INH", "KAN", "LZD", "MXF", "PZA", "RIF")
drug_list_val <- drug_list

lesion_level <- c("Plasma", "Lung", "Cellular lesion", "Caseum", "Cavity wall", "Caseous lesion")

COLOR_map <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
  "#CC6677", "#882255", "#AA4499", "#DDDDDD", "#332288", "#88CCEE",
  "#44AA99", "#117733"
)

paths <- list(
  data = here("data", "derived_datasets"),
  results = here("results")
)

dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$source_data, recursive = TRUE, showWarnings = FALSE)

# PowerPoint export: one deck, one slide per figure (graph2ppt from package export)
ppt_path <- file.path(paths$results, "preclinical_figures.pptx")
.ppt_append <- FALSE

#' Append a figure to `preclinical_figures.pptx` via export::graph2ppt().
#' Use `x` for ggplot objects; use `fun` for grid/base drawing (e.g. grid.arrange).
ppt_add_slide <- function(x = NULL, fun = NULL, width, height, ...) {
  if (!requireNamespace("export", quietly = TRUE)) {
    return(invisible(NULL))
  }
  if (is.null(x) && is.null(fun)) {
    stop("ppt_add_slide: provide either x (ggplot) or fun (drawing function).")
  }
  if (!is.null(x) && !is.null(fun)) {
    stop("ppt_add_slide: pass only one of x or fun.")
  }
  export::graph2ppt(
    x = x,
    fun = fun,
    file = ppt_path,
    width = width,
    height = height,
    append = .ppt_append,
    ...
  )
  .ppt_append <<- TRUE
  invisible(ppt_path)
}
