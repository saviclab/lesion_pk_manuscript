# Expects globals from 01_setup.R : %notin%, paths
# Expects sim_df already in environment from 06_simulation_coverage_plots.R
# - Panel A - known-outcome regimens (green = successful, grey = unsuccessful)
# - Panel B - includes predicted regimens (teal)

suppressPackageStartupMessages({
  library(cowplot)
})

# ---------------------------------------------------------------------------
# Step 1: Build per-drug coverage hours (Caseum / MRT)
# ---------------------------------------------------------------------------

coverage_all_drugs_REG_PLOTS <- sim_df %>%
  mutate(II_num = ifelse(II == "BID", 12, 24)) %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II_num, sep = "_")) %>%
  filter(LESIONNAME %in% c("Caseum")) %>%
  filter(LABEL == "MRT_mgL") %>%
  mutate(coverage_above_MRT = ifelse(CLESION > TARGET, 1, 0)) %>%
  group_by(ROW_ID) %>%
  summarise(hours_total_per_ROW_ID = sum(coverage_above_MRT),
            n_timepoints = n(),
            .groups = "drop") %>%
  mutate(drug_fraction = hours_total_per_ROW_ID / n_timepoints)

message("Max timepoints per drug: ", max(coverage_all_drugs_REG_PLOTS$n_timepoints))

message("coverage_all_drugs_REG_PLOTS ROW_IDs: ",
        paste(sort(coverage_all_drugs_REG_PLOTS$ROW_ID), collapse = ", "))

# ---------------------------------------------------------------------------
# Step 2: Known clinical regimen compositions
# ---------------------------------------------------------------------------

BPaMZ <- c("BDQ_400_mg_24",  "PMD_200_mg_24",  "MXF_400_mg_24", "PZA_1500_mg_24")
PaMZ  <- c("PMD_200_mg_24",  "MXF_400_mg_24",  "PZA_1500_mg_24")
HPZM  <- c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "MXF_400_mg_24")
HPZE  <- c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24")
HRZE  <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24")
HPZEC <- c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24", "CFZ_100_mg_24")
BLZ   <- c("BDQ_400_mg_24",  "LZD_600_mg_24",  "PZA_1500_mg_24")
BPaL  <- c("BDQ_400_mg_24",  "PMD_200_mg_24",  "LZD_600_mg_24")
BPaLM <- c("BDQ_400_mg_24",  "PMD_200_mg_24",  "LZD_600_mg_24", "MXF_400_mg_24")
RZME  <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "EMB_1200_mg_24")
HRZM  <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "INH_300_mg_24")
R15HZ <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "INH_300_mg_24")
R20HZ <- c("RIF_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24")
R40HZ <- c("RIF_1800_mg_24", "PZA_1500_mg_24", "INH_300_mg_24")
RZM   <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24")
HRZG  <- c("RIF_600_mg_24",  "PZA_1500_mg_24", "GTX_400_mg_24", "INH_300_mg_24")

# ---------------------------------------------------------------------------
# Step 3: regimen_df (16 known-outcome regimens)
# ---------------------------------------------------------------------------

make_reg_row <- function(row_ids, label) {
  coverage_all_drugs_REG_PLOTS %>%
    filter(ROW_ID %in% row_ids) %>%
    summarise(hours_total = sum(hours_total_per_ROW_ID),
              no_of_drugs = sum(drug_fraction)) %>%
    mutate(Regimen = label)
}

regimen_df <- bind_rows(
  make_reg_row(BPaMZ, "BPaMZ | Simplici-TB"),
  make_reg_row(HRZG,  "HRZG | OFLOTUB"),
  make_reg_row(PaMZ,  "PaMZ | STAND"),
  make_reg_row(HPZM,  "PMZ | S31/A5349, with H"),
  make_reg_row(HPZE,  "HPZE | S31/A5349"),
  make_reg_row(HRZE,  "HRZE | Standard-of-care"),
  make_reg_row(HPZEC, "HPZEC | CLO-FAST"),
  make_reg_row(BLZ,   "BLZ"),
  make_reg_row(BPaL,  "BPaL | Nix/ZeNix"),
  make_reg_row(BPaLM, "BPaLM | TB-PRACTECAL"),
  make_reg_row(RZME,  "RZME | REMoxTB"),
  make_reg_row(HRZM,  "HRZM | REMoxTB"),
  make_reg_row(R15HZ, "R15HZ | HIRIF/RIFASHORT"),
  make_reg_row(R20HZ, "R20HZ | HIRIF/RIFASHORT"),
  coverage_all_drugs_REG_PLOTS %>%
    filter(ROW_ID %in% R40HZ) %>%
    summarise(hours_total = sum(hours_total_per_ROW_ID) + 8,
              no_of_drugs = sum(drug_fraction) + 8 / max(coverage_all_drugs_REG_PLOTS$n_timepoints)) %>%
    mutate(Regimen = "R40HZ | HIRIF/RIFASHORT"),
  make_reg_row(RZM,   "RZM | REMoxTB")
)

lesion_ranked_regimens <- regimen_df %>%
  rename(Predictor = Regimen) %>%
  mutate(Predictor = ifelse(Predictor == "PMZ | S31/A5349, with H",
                            "HPZM | S31/A5349", Predictor)) %>%
  mutate(Outcome = ifelse(hours_total >= 48 & Predictor != "PaMZ | STAND", 1, 0))

clinical_regimen_outcome <- lesion_ranked_regimens %>% select(Predictor, Outcome)

message("Regimen coverage check (hours_total | no_of_drugs):")
for (i in seq_len(nrow(lesion_ranked_regimens))) {
  r <- lesion_ranked_regimens[i, ]
  message(sprintf("  %-30s hours=%3d  no_of_drugs=%.3f  Outcome=%d",
                  r$Predictor, r$hours_total, r$no_of_drugs, r$Outcome))
}

# ---------------------------------------------------------------------------
# Step 4: Full regimens list (known + novel + predicted)
#   Drug name updates: OPC-167832 -> QBS, GSK-656 -> GFB
# ---------------------------------------------------------------------------

regimens <- list(
  BPaMZ = c("BDQ_400_mg_24",  "PMD_200_mg_24",  "MXF_400_mg_24", "PZA_1500_mg_24"),
  PaMZ  = c("PMD_200_mg_24",  "MXF_400_mg_24",  "PZA_1500_mg_24"),
  HPZM  = c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "MXF_400_mg_24"),
  HPZE  = c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24"),
  HRZE  = c("RIF_600_mg_24",  "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24"),
  HPZEC = c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24", "CFZ_100_mg_24"),
  BLZ   = c("BDQ_400_mg_24",  "LZD_600_mg_24",  "PZA_1500_mg_24"),
  BPaL  = c("BDQ_400_mg_24",  "PMD_200_mg_24",  "LZD_600_mg_24"),
  BPaLM = c("BDQ_400_mg_24",  "PMD_200_mg_24",  "LZD_600_mg_24", "MXF_400_mg_24"),
  RZME  = c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "EMB_1200_mg_24"),
  HRZM  = c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "INH_300_mg_24"),
  HRZG  = c("RIF_600_mg_24",  "PZA_1500_mg_24", "GTX_400_mg_24", "INH_300_mg_24"),

  R15HZ     = c("RIF_600_mg_24",  "PZA_1500_mg_24", "INH_300_mg_24"),
  R20HZ     = c("RIF_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24"),
  R40HZ     = c("RIF_1800_mg_24", "PZA_1500_mg_24", "INH_300_mg_24"),
  RZM       = c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24"),

  R1200HZE  = c("RIF_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24"),
  R1800HZE  = c("RIF_1800_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24"),

  BMZ       = c("BDQ_400_mg_24", "MXF_400_mg_24", "PZA_1500_mg_24"),
  BPaUG8    = c("BDQ_400_mg_24", "PMD_200_mg_24", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BPaUG2    = c("BDQ_400_mg_24", "PMD_200_mg_24", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BOUG2     = c("BDQ_400_mg_24", "QBS_30_mg_24",  "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BPaOU     = c("BDQ_400_mg_24", "PMD_200_mg_24", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BDOG2     = c("BDQ_400_mg_24", "DLM_100_mg_12", "QBS_30_mg_24"),
  BDUG2     = c("BDQ_400_mg_24", "DLM_100_mg_12", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BPaOG2    = c("BDQ_400_mg_24", "PMD_200_mg_24", "QBS_30_mg_24"),
  BDOU      = c("BDQ_400_mg_24", "DLM_100_mg_12", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  DBO       = c("BDQ_400_mg_24", "DLM_100_mg_12", "QBS_30_mg_24"),
  BDUG8     = c("BDQ_400_mg_24", "DLM_100_mg_12", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  PaOUG2    = c("PMD_200_mg_24", "QBS_30_mg_24",  "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  DOUG2     = c("DLM_100_mg_12", "QBS_30_mg_24",  "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BDOG8     = c("BDQ_400_mg_24", "DLM_100_mg_12", "QBS_30_mg_24"),

  S8_25_PaL  = c("TBAJ-876_25_mg_24",  "PMD_200_mg_24", "LZD_600_mg_24"),
  S8_50_PaL  = c("TBAJ-876_50_mg_24",  "PMD_200_mg_24", "LZD_600_mg_24"),
  S8_100_PaL = c("TBAJ-876_100_mg_24", "PMD_200_mg_24", "LZD_600_mg_24"),

  R40MHZ        = c("RIF_1800_mg_24",     "MXF_400_mg_24", "PZA_1500_mg_24", "INH_300_mg_24"),
  S8_100_PaOU   = c("TBAJ-876_100_mg_24", "PMD_200_mg_24", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S5_100_PaOU   = c("TBAJ-587_100_mg_24", "PMD_200_mg_24", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S8_100_DOU    = c("TBAJ-876_100_mg_24", "DLM_100_mg_12", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S5_100_DOU    = c("TBAJ-587_100_mg_24", "DLM_100_mg_12", "QBS_30_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),

  S8_100_PaG6U  = c("TBAJ-587_100_mg_24", "PMD_200_mg_24", "GFB_20_mg_24",   "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S8_100_DG6U   = c("TBAJ-587_100_mg_24", "DLM_100_mg_12", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BPaG6U        = c("BDQ_400_mg_24",      "PMD_200_mg_24", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  BPaG6G2       = c("BDQ_400_mg_24",      "PMD_200_mg_24", "GFB_20_mg_24"),

  BMZRb  = c("BDQ_200_mg_24", "MXF_400_mg_24", "PZA_1500_mg_24", "RBT_300_mg_24"),
  BMZD   = c("BDQ_200_mg_24", "MXF_400_mg_24", "PZA_1500_mg_24", "DLM_100_mg_12"),

  S8MZRb  = c("TBAJ-876_100_mg_24", "PZA_1500_mg_24", "MXF_400_mg_24", "RBT_300_mg_24"),
  S8ZMD   = c("TBAJ-876_100_mg_24", "PZA_1500_mg_24", "MXF_400_mg_24", "DLM_100_mg_12"),
  S8ZMU   = c("TBAJ-876_100_mg_24", "PZA_1500_mg_24", "MXF_400_mg_24", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S8PaMU  = c("TBAJ-876_100_mg_24", "PMD_200_mg_24",  "MXF_400_mg_24", "SZD_1200_mg_24", "SZD-M1_1200_mg_24"),
  S8PaMRb = c("TBAJ-876_100_mg_24", "PMD_200_mg_24",  "MXF_400_mg_24", "RBT_300_mg_24"),
  S8PaMO  = c("TBAJ-876_100_mg_24", "PMD_200_mg_24",  "MXF_400_mg_24", "QBS_30_mg_24"),

  BDM    = c("BDQ_200_mg_24", "DLM_100_mg_12", "MXF_400_mg_24"),
  BDMG6  = c("BDQ_200_mg_24", "DLM_100_mg_12", "MXF_400_mg_24",  "GFB_20_mg_24"),
  BDZG6  = c("BDQ_200_mg_24", "DLM_100_mg_12", "PZA_1500_mg_24", "GFB_20_mg_24"),
  BDLG6  = c("BDQ_200_mg_24", "DLM_100_mg_12", "LZD_600_mg_24",  "GFB_20_mg_24"),
  BPaMG6 = c("BDQ_400_mg_24", "PMD_200_mg_24", "MXF_400_mg_24",  "GFB_20_mg_24"),
  BDMT   = c("BDQ_200_mg_24", "DLM_100_mg_12", "MXF_400_mg_24",  "BTZ-043_1250_mg_24"),
  BDZT   = c("BDQ_200_mg_24", "DLM_100_mg_12", "PZA_1500_mg_24", "BTZ-043_1250_mg_24"),
  BDLT   = c("BDQ_200_mg_24", "DLM_100_mg_12", "LZD_600_mg_24",  "BTZ-043_1250_mg_24"),
  BPaMT  = c("BDQ_400_mg_24", "PMD_200_mg_24", "MXF_400_mg_24",  "BTZ-043_1250_mg_24"),
  BMZT   = c("BDQ_200_mg_24", "MXF_400_mg_24", "PZA_1500_mg_24", "BTZ-043_1250_mg_24"),
  BDG6T  = c("BDQ_200_mg_24", "DLM_100_mg_12", "GFB_20_mg_24",   "BTZ-043_1250_mg_24")
)

# ---------------------------------------------------------------------------
# Step 4: regimen_coverage_all (loop over all regimens)
# ---------------------------------------------------------------------------

regimen_coverage_all <- NULL
for (i in seq_along(regimens)) {
  regimen_coverage_all <- bind_rows(
    regimen_coverage_all,
    coverage_all_drugs_REG_PLOTS %>%
      filter(ROW_ID %in% regimens[[i]]) %>%
      summarise(hours_total = sum(hours_total_per_ROW_ID),
                no_of_drugs = sum(drug_fraction)) %>%
      mutate(Regimen = names(regimens)[i])
  )
}

# ---------------------------------------------------------------------------
# Step 5: Panel A — known-outcome regimens only
# ---------------------------------------------------------------------------

p_fig6a <- lesion_ranked_regimens %>%
  mutate(Predictor = ifelse(Predictor == "R20HZ | HIRIF/RIFASHORT",
                            "R1200HZ | HIRIF/RIFASHORT", Predictor)) %>%
  mutate(Predictor = ifelse(Predictor == "R40HZ | HIRIF/RIFASHORT",
                            "R1800HZ | HIRIF/RIFASHORT", Predictor)) %>%
  filter(Predictor != "R15HZ | HIRIF/RIFASHORT") %>%
  filter(Predictor != "BLZ") %>%
  filter(Predictor != "RZM | REMoxTB") %>%
  mutate(no_of_drugs_plot = ifelse(no_of_drugs > 0, pmax(no_of_drugs, 0.05), 0)) %>%
  arrange(Outcome, no_of_drugs) %>%
  mutate(Predictor = fct_inorder(Predictor)) %>%
  ggplot(aes(x = no_of_drugs_plot, y = Predictor, fill = factor(Outcome))) +
  geom_col() +
  scale_fill_manual(values = rev(c("#87AE73", "grey"))) +
  guides(fill = "none") +
  ylab("") +
  xlab("Number of drugs in the regimen\nthat meet the PK-PD target in caseum") +
  theme_bw() +
  theme(axis.text.y  = element_text(size = 8),
        axis.title.x = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Step 6: Panel B — with predicted regimens
# ---------------------------------------------------------------------------

unite4tb <- c("BDM", "BDMG6", "BDZG6", "BDLG6", "BPaMG6",
              "BDMT", "BDZT", "BDLT", "BPaMT", "BMZT", "BDG6T")
nc009    <- c("S8_25_PaL", "S8_50_PaL", "S8_100_PaL")
crushtb  <- c("BMZRb", "BMZD")
pantb    <- c("BPaOU", "BDOU")

lesion_ranked_regimens_new_reg <- regimen_coverage_all %>%
  rename(regimen_name = Regimen) %>%
  mutate(regimen_name = ifelse(regimen_name == "R40MHZ", "R1800MHZ", regimen_name)) %>%
  filter(regimen_name %notin% c("BPaMZ", "PaMZ", "HPZM", "HPZE", "HRZE",
                                "HPZEC", "BPaL", "BPaLM", "BPaUG8", "R15HZ", "R20HZ", "R40HZ",
                                "R1200HZE", "R1800HZE", "RZME", "RZM", "HRZM")) %>%
  filter(!grepl("G2", regimen_name)) %>%
  filter(!grepl("G8", regimen_name)) %>%
  mutate(Predictor = regimen_name) %>%
  mutate(Outcome = -1) %>%
  bind_rows(lesion_ranked_regimens) %>%
  filter(Outcome >= 0 | Predictor %in% c(unite4tb, crushtb, nc009, pantb)) %>%
  mutate(Predictor = gsub("S8", "TBAJ876", Predictor)) %>%
  mutate(Predictor = gsub("S5_100_", "TBAJ587_", Predictor)) %>%
  mutate(Predictor = ifelse(Predictor == "R20HZ | HIRIF/RIFASHORT",
                            "R1200HZ | HIRIF/RIFASHORT", Predictor)) %>%
  mutate(Predictor = ifelse(Predictor == "R40HZ | HIRIF/RIFASHORT",
                            "R1800HZ | HIRIF/RIFASHORT", Predictor)) %>%
  filter(Predictor != "R15HZ | HIRIF/RIFASHORT") %>%
  filter(Predictor != "RZM | REMoxTB") %>%
  filter(Predictor != "BLZ")

fig6b_data <- lesion_ranked_regimens_new_reg %>%
  arrange(hours_total, Outcome, desc(Predictor)) %>%
  group_by(hours_total, Outcome) %>%
  mutate(Predictor = ifelse(Predictor == "BDOU",  "BDQU",  Predictor)) %>%
  mutate(Predictor = ifelse(Predictor == "BPaOU", "BPaQU", Predictor)) %>%
  ungroup() %>%
  mutate(tbaj_boost = case_when(
    Predictor == "TBAJ876_100_PaL" ~ 0.03,
    Predictor == "TBAJ876_50_PaL"  ~ 0.02,
    Predictor == "TBAJ876_25_PaL"  ~ 0.01,
    TRUE ~ 0
  )) %>%
  mutate(no_of_drugs_plot = ifelse(no_of_drugs > 0, pmax(no_of_drugs, 0.05), 0)) %>%
  mutate(sort_val = no_of_drugs + Outcome * 0.001 + tbaj_boost)

message("\nPanel B regimen order (bottom to top on plot):")
fig6b_data %>%
  arrange(sort_val) %>%
  mutate(row = row_number()) %>%
  {for (i in seq_len(nrow(.))) {
    r <- .[i, ]
    message(sprintf("  %2d  %-35s  no_of_drugs=%.3f  Outcome=%2d  sort=%.3f",
                    r$row, r$Predictor, r$no_of_drugs, r$Outcome, r$sort_val))
  }}

p_fig6b <- fig6b_data %>%
  mutate(Predictor = fct_reorder(Predictor, sort_val)) %>%
  ggplot(aes(x = no_of_drugs_plot, y = Predictor, fill = factor(Outcome))) +
  geom_col() +
  scale_fill_manual(values = c("#008080", "grey", "#87AE73")) +
  guides(fill = "none") +
  ylab("") +
  xlab("Number of drugs in the regimen\nthat meet the PK-PD target in caseum") +
  theme_bw() +
  annotate(geom = "text", label = "'Successful' regimen",   color = "#87AE73", x = 3, y = "BDLT",  hjust = 1) +
  annotate(geom = "text", label = "'Unsuccessful' regimen", color = "grey",    x = 3, y = "BDLG6", hjust = 1) +
  annotate(geom = "text", label = "Predicted regimen",      color = "#008080", x = 3, y = "BDG6T", hjust = 1) +
  theme(axis.text     = element_text(size = 8),
        axis.title.x  = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Step 7: Combine into Figure 6
# ---------------------------------------------------------------------------

p_figure6 <- cowplot::plot_grid(
  p_fig6a, p_fig6b,
  ncol = 1,
  labels = c("A.", "B."),
  rel_heights = c(2.5, 4.5)
)

fig6_path <- file.path(paths$results, "Fig6_regimen_coverage.png")
ggsave(fig6_path, p_figure6, width = 7, height = 7, dpi = 300)
message("Saved: ", fig6_path)

sim_ppt_add_slide(x = p_figure6, width = 7, height = 7)
message("Figure 6 appended to PPTX")
