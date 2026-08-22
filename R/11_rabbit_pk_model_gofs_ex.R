# Example code for outputting Goodness-of-fit plots

suppressPackageStartupMessages({
  library(vpc)
  library(xpose)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(reshape)
  library(export)
})

# drugs <- list of drug names
# WD <- main working directory
# wd_values <- list of file paths to NONMEM runs
# runno_values <- list of NONMEM final models

# ---------------------------------------------------------------------------
# Plasma
# ---------------------------------------------------------------------------

p_plasma_plot_list <- list() # create an empty list to store the plots

for (i in seq_along(drugs)) { # to access corresponding elements of each list for each loop
  drug <- drugs[i]
  wd <- wd_values[i]
  runno <- runno_values[i]
  
  setwd(paste0(WD, wd))
  xpdb <- xpose_data(runno = runno)
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_plasma <- dv_vs_ipred(xpdb, 
                          guide = T, 
                          type = "pl",
                          caption = NULL) +
    theme_bw(base_size = 16) +
    labs(title = drug, subtitle = "Plasma",
         x = "Individual predictions",
         y = "Observations")
  
  p_plasma_plot_list[[i]] <- p_plasma # store each plot in the list
}

for (i in seq_along(p_plasma_plot_list)) { # display each plot
  print(p_plasma_plot_list[[i]])
}

# ---------------------------------------------------------------------------
# Uninvolved Lung
# ---------------------------------------------------------------------------

p_lung_plot_list <- list()

for (i in seq_along(drugs)) {
  drug <- drugs[i]
  wd <- wd_values[i]
  runno_lung <- runnos_lung[i]
  
  setwd(paste0(WD, wd))
  xpdb <- xpose_data(runno = runno_lung)
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_lung <- dv_vs_ipred(xpdb, 
                        guide = T, 
                        type = "pl",
                        caption = NULL) +
    theme_bw(base_size = 16) +
    labs(title = "", subtitle = "Lung",
         x = "Individual predictions",
         y = "Observations")
  p_lung 
  
  p_lung_plot_list[[i]] <- p_lung
}

for (i in seq_along(p_lung_plot_list)) {
  print(p_lung_plot_list[[i]])
}

# ---------------------------------------------------------------------------
# Cellular Lesion
# ---------------------------------------------------------------------------

p_cell_plot_list <- list()

for (i in seq_along(drugs)) {
  drug <- drugs[i]
  wd <- wd_values[i]
  runno_cell <- runnos_cell[i]
  
  setwd(paste0(WD, wd))
  xpdb <- xpose_data(runno = runno_cell)
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_cell <- dv_vs_ipred(xpdb, 
                        guide = T, 
                        type = "pl",
                        caption = NULL) +
    theme_bw(base_size = 16) +
    labs(title = "", subtitle = "Cellular Lesion",
         x = "Individual predictions",
         y = "Observations")
  
  p_cell_plot_list[[i]] <- p_cell
}

for (i in seq_along(p_cell_plot_list)) {
  print(p_cell_plot_list[[i]])
}

# ---------------------------------------------------------------------------
# Caseum
# ---------------------------------------------------------------------------

p_case_plot_list <- list()

for (i in seq_along(drugs)) {
  drug <- drugs[i]
  wd <- wd_values[i]
  runno_case <- runnos_case[i]
  
  setwd(paste0(WD, wd))
  xpdb <- xpose_data(runno = runno_case)
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_case <- dv_vs_ipred(xpdb, 
                        guide = T, 
                        type = "pl",
                        caption = NULL) +
    theme_bw(base_size = 16) +
    labs(title = "", subtitle = "Caseum",
         x = "Individual predictions",
         y = "Observations")
  
  p_case_plot_list[[i]] <- p_case
}

for (i in seq_along(p_case_plot_list)) {
  print(p_case_plot_list[[i]])
}

# ---------------------------------------------------------------------------
# Arrange plots and save
# ---------------------------------------------------------------------------

a <- grid.arrange(p_plasma_plot_list[[1]], p_lung_plot_list[[1]], p_cell_plot_list[[1]], p_case_plot_list[[1]],
                  p_plasma_plot_list[[2]], p_lung_plot_list[[2]], p_cell_plot_list[[2]], p_case_plot_list[[2]],
                  p_plasma_plot_list[[3]], p_lung_plot_list[[3]], p_cell_plot_list[[3]], p_case_plot_list[[3]],
                  p_plasma_plot_list[[4]], p_lung_plot_list[[4]], p_cell_plot_list[[4]], p_case_plot_list[[4]],
                  ncol = 4)

graph2ppt({grid.draw(a)}, file = paste0(WD, "drug1_drug2_drug3_drug4.pptx"), width = 15, height = 12)

