# Example code for outputting Visual Predictive Checks

suppressPackageStartupMessages({
  library(vpc)
  library(xpose)
  library(xpose4)
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
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_plasma <- vpc::vpc_vpc(psn_folder = paste0(WD, 
                                               wd,
                                               'run', 
                                               runno, 
                                               'vpc-idv'),
                           obs_cols = list(idv = "TAD"),
                           sim_cols = list(idv = "TAD"),
                           bins = "pretty",
                           pred_corr = F,
                           # logy = T,
                           show = list(obs_dv = T, obs_ci = F, pi = T, pi_as_area = T)) +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank()) +
    scale_x_continuous(breaks = seq(0,120,12)) +
    # scale_y_log10() +
    labs(title = drug, subtitle = "Plasma",
         x = "Time after dose (h)",
         y = "Concentration (mg/L)")
  
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
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_lung <- vpc::vpc_vpc(psn_folder = paste0(WD, 
                                             wd,
                                             'Lung/',
                                             'vpc_',
                                             'run', 
                                             runno_lung),
                         obs_cols = list(idv = "TAD"),
                         sim_cols = list(idv = "TAD"),
                         bins = "pretty",
                         pred_corr = F,
                         show = list(obs_dv = T, obs_ci = F, pi = T, pi_as_area = T)) +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank()) +
    scale_x_continuous(breaks = seq(0,120,12)) +
    labs(title = "", subtitle = "Lung",
         x = "Time after dose (h)",
         y = "Concentration (mg/kg)")
  
  p_lung_plot_list[[i]] <- p_lung
}

for (i in seq_along(p_lung_plot_list)) {
  print(p_lung_plot_list[[i]])
}

# ---------------------------------------------------------------------------
# Cellular lesion
# ---------------------------------------------------------------------------

p_cell_plot_list <- list()

for (i in seq_along(drugs)) {
  drug <- drugs[i]
  wd <- wd_values[i]
  runno_cell <- runnos_cell[i]
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_cell <- vpc::vpc_vpc(psn_folder = paste0(WD, 
                                             wd,
                                             'Cell/',
                                             'vpc_',
                                             'run', 
                                             runno_cell),
                         obs_cols = list(idv = "TAD"),
                         sim_cols = list(idv = "TAD"),
                         bins = "pretty",
                         pred_corr = F,
                         show = list(obs_dv = T, obs_ci = F, pi = T, pi_as_area = T)) +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank()) +
    scale_x_continuous(breaks = seq(0,120,12)) +
    labs(title = "", subtitle = "Cellular Lesion",
         x = "Time after dose (h)",
         y = "Concentration (mg/kg)")
  
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
  
  options(repr.plot.width = 6, repr.plot.height = 4)
  p_case <- vpc::vpc_vpc(psn_folder = paste0(WD, 
                                             wd,
                                             'Caseum/',
                                             'vpc_',
                                             'run', 
                                             runno_case),
                         obs_cols = list(idv = "TAD"),
                         sim_cols = list(idv = "TAD"),
                         bins = "pretty",
                         pred_corr = F,
                         show = list(obs_dv = T, obs_ci = F, pi = T, pi_as_area = T)) +
    theme_bw(base_size = 16) +
    theme(panel.grid = element_blank()) +
    scale_x_continuous(breaks = seq(0,120,12)) +
    labs(title = "", subtitle = "Caseum",
         x = "Time after dose (h)",
         y = "Concentration (mg/kg)")
  
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
