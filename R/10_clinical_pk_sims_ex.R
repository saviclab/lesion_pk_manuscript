# Example code for clinical simulations in PKPDsim

suppressPackageStartupMessages({
  library(tidyverse)
  library(xpose)
  library(xpose4)
  library(ggpubr)
  library(gridExtra)
  library(grid)
  library(readxl)
  library(ggsci)
  library(PKPDsim)
})

# ---------------------------------------------------------------------------
# Example model structure: 1 cmt, 1 transit
# ---------------------------------------------------------------------------

ex_model <- new_ode_model(code = "dAdt[1] = -KTR*A[1] 
                                  dAdt[2] = KA*A[3] - CL/V*A[2] 
                                  dAdt[3] = KTR*A[1] - KA*A[3]
                                  CP = A[2]/V
                                  dAdt[4] = KPL * (PC*A[2]/V - A[4])
                                  CLESION = A[4]
                                  dAdt[5] = A[2]/V
                                  dAdt[6] = A[4]
                                  AUC = A[5]
                                  AUCL = A[6]
                                 ",
                           obs = list(cmt = 2, scale = "V"),
                           dose = list(cmt = 1),
                           declare_variables = c("CP","CLESION", "AUC", "AUCL"))

ex_model

# ---------------------------------------------------------------------------
# Define regimen
# ---------------------------------------------------------------------------

ex_regimen <- new_regimen(amt = 100,
                          interval = 24, 
                          n = 28) 

# ---------------------------------------------------------------------------
# Define parameters
# ---------------------------------------------------------------------------

KTR = a
KA = b
CL = c
V = d

# human
ex_PC_vector =  c(e, f, g)
ex_KPL_vector = c(h, i, j)

KPL_par_vector = c("kpl_lung", "kpl_cell", "kpl_caseum")
PC_par_vector = c("PC_lung", "PC_cell", "PC_caseum")
lesionname_vector = c("Lung", "Cellular Lesion", "Caseum")

# ---------------------------------------------------------------------------
# Run simulation
# ---------------------------------------------------------------------------

ex_simulation_collect <- data.frame()

for (j in 1:length(ex_KPL_vector)) {
  
  ex_parameters <- list("KTR" = KTR,
                        "KA" = KA,
                        "CL" = CL,
                        "V" = V,
                        "KPL" = ex_KPL_vector[j],
                        "PC" = ex_PC_vector[j]
  )
  
  ex_simulation <- sim(ode = ex_model,
                       regimen = ex_regimen,
                       parameters = ex_parameters,
                       only_obs = T,
                       output_include = list(variables=T),
                       t_obs = seq(0, 14*24)
  ) %>% 
    mutate(DRUG = "Drug Name",
           DOSE = 100,
           II = "QD",
           TIME = t,
           LESIONNAME = lesionname_vector[j]
    ) %>% 
    select(DRUG, DOSE, II, LESIONNAME, TIME, CP, CLESION, AUC, AUCL)
  names(ex_simulation)
  ex_simulation_collect <- rbind(ex_simulation_collect, ex_simulation)
}
