# set up
library(tidyverse)
library(xpose)
library(ggpubr)
library(gridExtra)
library(grid)
library(readxl)
library(ggsci)
library(PKPDsim)
library(pracma) 

# model
my_model_oral <- new_ode_model(code = "tad = t - t_prv_dose
                                       dAdt[1] = -KA*A[1] 
                                       dAdt[2] = KA*A[1] - CL/V*A[2]
                                       dAdt[3] = A[2]/V
                                       CP = A[2]/V
                                       dAdt[4] = KPL * (PC*A[2]/V - A[4])
                                       CLESION = A[4]
                                       dAdt[5] = A[4]
                                       AUC = A[3]
                                       AUCL = A[5]",
                               declare_variables = c("tad", "AUC", "AUCL", "CLESION", "CP"),
                               dose = list(cmt = 1),
                               obs = list(cmt = 2, scale = "V")
                              )
my_model_oral

# parameters
my_parameters_oral <- list(CL = 4.97,
                           V = 121.56, 
                           KA = 0.29,
                           KPL = 0.19,
                           PC = 2.08)

# regimen
my_regimen <- new_regimen(amt = 30,
                          interval =24, 
                          n = 24) 

# simulate
my_simulated_df_oral <- sim(ode = my_model_oral,
                            regimen = my_regimen,
                            parameters = my_parameters_oral,
                            only_obs = T,
                            output_include = list(variables=T),
                            t_obs = seq(0,12*24))

steady_state_time <- tail(my_simulated_df_oral$t, 24)
steady_state_concentration <- tail(my_simulated_df_oral$CP, 24)
auc_0_24 <- trapz(steady_state_time, steady_state_concentration)
CMAX <- max(steady_state_concentration)
TMAX <- steady_state_time[which.max(steady_state_concentration)]

df_to_plot <- my_simulated_df_oral %>% 
  filter(t >= max(t)-24) %>% 
  mutate(TAD = t - (max(t)-24)) %>% 
  filter(CLESION != 0)
# df_to_plot$CLESION>0.85

curve_p <- ggplot(df_to_plot, mapping = aes(x = TAD, y = CP)) +
  geom_line(linewidth = 1) +
  labs(y = "Concentration (mg/L)", 
       x = "Time (hours)", 
       title = paste("TBD-11", "30", "mg", sep=" ")) +
  scale_x_continuous(breaks = seq(0,24,4)) +
  theme_classic()
curve_p

curve_p + 
  geom_hline(aes(yintercept = 0.85, color =  'casMBC90'))

df_to_plot$coverage_above_MRT = ifelse(df_to_plot$CLESION>0.85, 1, 0)

# need to fix this ggplot (if no MRT==0 in dataset, then results are incorrect)
ggplot(df_to_plot, mapping = aes(x = TAD, y = comp, fill = factor(coverage_above_MRT))) +
  geom_point(shape = 22, size = 4, stroke = 0.2) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0,24,4)) +
  # facet_wrap(~LESIONNAME) +
  labs(y = "Drug",
       x = "Time after dose (hours)") +
  scale_fill_manual(values = c("navy","white"), 
                    name = "",
                    labels = c("above target","below target"))
