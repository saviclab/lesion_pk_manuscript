# Self-contained: only requires here() paths and 3 CSV data files
# Expects sim_ppt_add_slide() to be in the global environment (from 06)

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

blue   <- "#0072B2"
orange <- "#D55E00"

# ---------------------------------------------------------------------------
# Load datasets
# ---------------------------------------------------------------------------

clin_data <- here("data", "derived_datasets", "clin_outcome_datasets", "final_data")

DF_TBRS31 <- read.csv(file.path(clin_data, "TBR_S31_18m.csv")) %>%
  filter(EVID == 0) %>%
  mutate(DV = as.numeric(DV),
         LUNGGRADING = -99)

DF_RIFASHORT <- read.csv(file.path(clin_data, "RIFASHORT_18m_wLUNGGRADING.csv")) %>%
  filter(EVID == 0) %>%
  mutate(DV = as.numeric(DV))

stopifnot(identical(colnames(DF_TBRS31), colnames(DF_RIFASHORT)))

DF_TBR <- DF_TBRS31 %>% filter(STUDYID != 2000)
DF_S31 <- DF_TBRS31 %>% filter(STUDYID == 2000)

# ---------------------------------------------------------------------------
# Risk stratification
# ---------------------------------------------------------------------------

# TBReflect studies (cavity + smear)
DF_TBR <- DF_TBR %>%
  filter(CAVITY != -99, SMEAR != -99) %>%
  mutate(GROUP = ifelse(CAVITY == 1 & SMEAR >= 2, "Hard", "Moderate/easy"),
         LUNGGRADING1 = -99)

# Study 31 (Xpert + CXR extent)
DF_S31 <- DF_S31 %>%
  filter(XPERT != -99, CXREXTNT != -99) %>%
  mutate(GROUP = ifelse(XPERT < 18 & CXREXTNT >= 2, "Hard", "Moderate/easy"),
         LUNGGRADING1 = -99)

# RIFASHORT (Xpert + lung grading)
DF_RIFASHORT <- DF_RIFASHORT %>%
  mutate(LUNGGRADING  = as.factor(LUNGGRADING),
         LUNGGRADING1 = as.factor(ifelse(is.na(LUNGGRADING), 1,
                                         ifelse(LUNGGRADING == 4, 1, 0))))

DF_RIFASHORT_XPERT <- DF_RIFASHORT %>%
  filter(XPERT != -99) %>%
  mutate(GROUP = ifelse(XPERT < 19 & LUNGGRADING1 == 1, "Hard", "Moderate/easy"))

DF_COMBINED_XPERT <- rbind(rbind(DF_TBR, DF_S31), DF_RIFASHORT_XPERT)

# ---------------------------------------------------------------------------
# Map STUDYID + ARM -> regimen labels
# ---------------------------------------------------------------------------

study_reg_factors <- c(
  "REMoxTB: 2HRZE/4HR", "REMoxTB: 2HRZM/2HRM", "REMoxTB: 2MRZE/2MR",
  "RIFAQUIN: 2HRZE/4HR", "RIFAQUIN: 2MRZE/2P2M2", "RIFAQUIN: 2MRZE/4P1M1",
  "OFLOTUB: 2HRZE/4HR", "OFLOTUB: 2HRZG/2HRG",
  "TBRU: 2HRZE/4HR", "TBRU: 2HRZE/2HR",
  "S31: 2HRZE/4HR", "S31: 2HPZE/2HP", "S31: 2HPZM/2HPM",
  "RIFASHORT: 2HRZE/4HR", "RIFASHORT: 2HRZE/2HR (R1200mg)",
  "RIFASHORT: 2HRZE/2HR (R1800mg)", "Nix/ZeNix: BPaL"
)

DF_KM_XPERT <- DF_COMBINED_XPERT %>%
  mutate(Status = as.numeric(DV)) %>%
  mutate(STUDY_REG = case_when(
    STUDYID == 1021 & ARM == 1  ~ "REMoxTB: 2HRZE/4HR",
    STUDYID == 1021 & ARM == 2  ~ "REMoxTB: 2HRZM/2HRM",
    STUDYID == 1021 & ARM == 3  ~ "REMoxTB: 2MRZE/2MR",
    STUDYID == 1020 & ARM == 1  ~ "RIFAQUIN: 2HRZE/4HR",
    STUDYID == 1020 & ARM == 5  ~ "RIFAQUIN: 2MRZE/2P2M2",
    STUDYID == 1020 & ARM == 6  ~ "RIFAQUIN: 2MRZE/4P1M1",
    STUDYID == 1022 & ARM == 1  ~ "OFLOTUB: 2HRZE/4HR",
    STUDYID == 1022 & ARM == 8  ~ "OFLOTUB: 2HRZG/2HRG",
    STUDYID == 1001 & ARM == 1  ~ "TBRU: 2HRZE/4HR",
    STUDYID == 1001 & ARM == 9  ~ "TBRU: 2HRZE/2HR",
    STUDYID == 2000 & ARM == 1  ~ "S31: 2HRZE/4HR",
    STUDYID == 2000 & ARM == 10 ~ "S31: 2HPZE/2HP",
    STUDYID == 2000 & ARM == 11 ~ "S31: 2HPZM/2HPM",
    STUDYID == 3000 & ARM == 1  ~ "RIFASHORT: 2HRZE/4HR",
    STUDYID == 3000 & ARM == 12 ~ "RIFASHORT: 2HRZE/2HR (R1200mg)",
    STUDYID == 3000 & ARM == 13 ~ "RIFASHORT: 2HRZE/2HR (R1800mg)",
    .default = NA_character_
  )) %>%
  select(ID, TIME, Status, RISK_GROUP = GROUP, STUDY_REG) %>%
  mutate(CENS = ifelse(Status == 0, "Yes: Yes", "No: No")) %>%
  mutate(STUDY_REG  = factor(STUDY_REG, levels = study_reg_factors)) %>%
  mutate(RISK_GROUP = factor(RISK_GROUP, levels = c("Moderate/easy", "Hard")))

# ---------------------------------------------------------------------------
# BPaL (Nix/ZeNix) — cavity + smear stratification
# ---------------------------------------------------------------------------

all_mdr <- read.csv(file.path(clin_data, "All_validationMDR.csv"), header = TRUE)
all_mdr$DV <- as.numeric(as.character(all_mdr$DV))
all_mdr <- all_mdr %>%
  rename(TIMEM = TIME2) %>%
  mutate(RISK_MDR4 = ifelse(CAVITY < 2 & SMEAR < 2, "Low", "Moderate/High"))

DF_BPAL <- all_mdr %>%
  group_by(SUBJID) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(TIME = TIMEM) %>%
  rename(Status = DV) %>%
  mutate(RISK_GROUP = ifelse(RISK_MDR4 == "Moderate/High", "Hard", "Moderate/easy")) %>%
  mutate(STUDY_REG = "Nix/ZeNix: BPaL") %>%
  mutate(CENS = ifelse(Status == 0, "Yes: Yes", "No: No")) %>%
  select(all_of(names(DF_KM_XPERT)))

# ---------------------------------------------------------------------------
# Combine & filter — combine RIFASHORT arms
# ---------------------------------------------------------------------------

ordered_study_reg_factors <- c(
  "Nix/ZeNix: BPaL",
  "S31: 2HPZM/2HPM",
  "RIFASHORT: 2HRZE/2HR (R1200mg & R1800mg)",
  "REMoxTB: 2HRZM/2HRM",
  "OFLOTUB: 2HRZG/2HRG",
  "REMoxTB: 2MRZE/2MR",
  "S31: 2HPZE/2HP"
)

ALL_DF <- bind_rows(
  DF_KM_XPERT %>% mutate(ID = as.character(ID)),
  DF_BPAL     %>% mutate(ID = as.character(ID))
) %>%
  filter(!grepl("RIFAQUIN", STUDY_REG),
         !grepl("TBRU",     STUDY_REG),
         !grepl("2HRZE/4HR", STUDY_REG)) %>%
  mutate(STUDY_REG = ifelse(grepl("RIFASHORT", STUDY_REG),
                            "RIFASHORT: 2HRZE/2HR (R1200mg & R1800mg)",
                            as.character(STUDY_REG))) %>%
  mutate(RISK_GROUP = factor(RISK_GROUP, levels = c("Moderate/easy", "Hard"))) %>%
  mutate(STUDY_REG  = factor(STUDY_REG, levels = ordered_study_reg_factors))

# ---------------------------------------------------------------------------
# KM fit + faceted plot
# ---------------------------------------------------------------------------

ALL_fit <- survfit(Surv(TIME, Status) ~ RISK_GROUP, data = ALL_DF)

ALL_plot <- ggsurvplot_facet(
  ALL_fit,
  data          = ALL_DF,
  facet.by      = "STUDY_REG",
  short.panel.labs = TRUE,
  palette       = c(blue, orange),
  conf.int      = FALSE,
  censor        = FALSE,
  pval          = FALSE,
  risk.table    = FALSE,
  legend.labs   = c("Easy-to-treat", "Hard-to-treat"),
  break.time.by = 6
)

ALL_plot <- ALL_plot +
  coord_cartesian(xlim = c(0, 18), ylim = c(0.75, 1)) +
  labs(x = "Time since treatment initiation (months)",
       y = "Proportion without TB-related\nunfavorable outcomes") +
  facet_wrap(~STUDY_REG,
             labeller = labeller(STUDY_REG = label_wrap_gen(14), .multi_line = FALSE),
             nrow = 1, scales = "fixed") +
  theme_bw() +
  theme(panel.grid       = element_blank(),
        strip.background = element_rect(fill = "transparent"),
        legend.position  = "top",
        legend.title     = element_blank())

df_n_patients <- as.data.frame(table(ALL_DF$STUDY_REG, ALL_DF$RISK_GROUP)) %>%
  mutate(STUDY_REG = Var1)

df_pvals <- ALL_DF %>%
  group_by(STUDY_REG) %>%
  summarise(pval = survdiff(Surv(TIME, Status) ~ RISK_GROUP, data = pick(everything()))$pvalue,
            .groups = "drop") %>%
  mutate(pval_label = paste("p =", signif(pval, 2)))

p_figure5 <- ALL_plot +
  geom_text(df_n_patients %>% filter(Var2 == "Moderate/easy"),
            mapping = aes(x = 1, y = 0.81, label = paste("n =", Freq)),
            hjust = 0, color = blue) +
  geom_text(df_n_patients %>% filter(Var2 == "Hard"),
            mapping = aes(x = 1, y = 0.79, label = paste("n =", Freq)),
            hjust = 0, color = orange) +
  geom_text(df_pvals,
            mapping = aes(x = 1, y = 0.77, label = pval_label),
            hjust = 0, color = "black")

# ---------------------------------------------------------------------------
# Panel B: Per-drug PK-PD coverage bars (single faceted plot for alignment)
#   Easy-to-treat = Cellular lesion concentration over macIC90 (= MRT)
#   Hard-to-treat = Caseum concentration over casMBC50 (= MRT)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(cowplot))

fig5_regimens <- list(
  "Nix/ZeNix: BPaL"                          = c("BDQ_400_mg_24",  "PMD_200_mg_24",  "LZD_600_mg_24"),
  "S31: 2HPZM/2HPM"                          = c("RPT_1200_mg_24", "PZA_1500_mg_24", "MXF_400_mg_24", "INH_300_mg_24"),
  "RIFASHORT: 2HRZE/2HR (R1200mg & R1800mg)" = c("RIF_1800_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24"),
  "REMoxTB: 2HRZM/2HRM"                      = c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "INH_300_mg_24"),
  "OFLOTUB: 2HRZG/2HRG"                      = c("RIF_600_mg_24",  "PZA_1500_mg_24", "GTX_400_mg_24", "INH_300_mg_24"),
  "REMoxTB: 2MRZE/2MR"                       = c("RIF_600_mg_24",  "PZA_1500_mg_24", "MXF_400_mg_24", "EMB_1200_mg_24"),
  "S31: 2HPZE/2HP"                           = c("RPT_1200_mg_24", "PZA_1500_mg_24", "INH_300_mg_24", "EMB_1200_mg_24")
)

totals_for_levels_fig5 <- coverage_all_drugs_intermediate %>%
  mutate(NEW_ROW_ID = paste(DRUG, DOSE, "mg", II_num, sep = "_")) %>%
  filter(LESIONNAME == "Caseum") %>%
  select(DRUG, NEW_ROW_ID, coverage_above_MRT) %>%
  group_by(DRUG, NEW_ROW_ID) %>%
  summarise(count = sum(coverage_above_MRT), .groups = "drop")

coverage_easy_hard <- bind_rows(
  coverage_all_drugs_intermediate %>%
    filter(LESIONNAME == "Cellular lesion", coverage_above_MRT == 1) %>%
    mutate(LESIONNAME2 = "Easy-to-treat"),
  coverage_all_drugs_intermediate %>%
    filter(LESIONNAME == "Caseum", coverage_above_MRT == 1) %>%
    mutate(LESIONNAME2 = "Hard-to-treat")
) %>%
  mutate(NEW_ROW_ID = paste(DRUG, DOSE, "mg", II_num, sep = "_"),
         LESIONNAME2 = factor(LESIONNAME2, levels = c("Easy-to-treat", "Hard-to-treat")))

# Uniform number of rows per facet so coverage is bottom-anchored and aligned
# Across regimens (pad shorter regimens with blank placeholder rows on top).
max_drugs <- max(lengths(fig5_regimens))

# Per-regimen drug ordering (most caseum coverage -> bottom) + top blank pads
reg_key_info <- purrr::map(names(fig5_regimens), function(reg_name) {
  drug_ids   <- fig5_regimens[[reg_name]]
  drug_names <- sub("_.*", "", drug_ids)

  drug_counts <- totals_for_levels_fig5 %>%
    filter(NEW_ROW_ID %in% drug_ids) %>%
    select(DRUG, count)
  all_counts <- tibble(DRUG = drug_names) %>%
    left_join(drug_counts, by = "DRUG") %>%
    mutate(count = replace_na(count, 0))
  ordered_drugs <- all_counts %>% arrange(desc(count)) %>% pull(DRUG)

  real_keys <- paste0(ordered_drugs, "___", reg_name)
  n_pad     <- max_drugs - length(real_keys)
  pad_keys  <- if (n_pad > 0) paste0(strrep(" ", seq_len(n_pad)), "___", reg_name) else character(0)

  list(keys = c(real_keys, pad_keys)) # Ordered bottom -> top
})
names(reg_key_info) <- names(fig5_regimens)

coverage_fig5 <- purrr::map_dfr(names(fig5_regimens), function(reg_name) {
  drug_ids <- fig5_regimens[[reg_name]]

  reg_data <- coverage_easy_hard %>%
    filter(NEW_ROW_ID %in% drug_ids) %>%
    mutate(STUDY_REG = reg_name,
           DRUG_KEY  = paste0(DRUG, "___", reg_name))

  # Scaffold every row slot (real drugs + blank pads) x both lesion panels
  scaffold <- expand_grid(
    DRUG_KEY    = reg_key_info[[reg_name]]$keys,
    LESIONNAME2 = factor(c("Easy-to-treat", "Hard-to-treat"),
                         levels = c("Easy-to-treat", "Hard-to-treat"))
  ) %>% mutate(STUDY_REG = reg_name, TAD = NA_real_)

  bind_rows(reg_data, scaffold)
})

drug_key_order <- purrr::map(ordered_study_reg_factors, function(reg_name) {
  reg_key_info[[reg_name]]$keys
}) %>% unlist()

coverage_fig5$DRUG_KEY <- factor(coverage_fig5$DRUG_KEY, levels = drug_key_order)

facet_levels_B <- as.vector(sapply(ordered_study_reg_factors, function(reg) {
  c(paste0("E|", reg), paste0("H|", reg))
}))

coverage_fig5 <- coverage_fig5 %>%
  mutate(facet_var = paste0(ifelse(LESIONNAME2 == "Easy-to-treat", "E|", "H|"),
                            STUDY_REG),
         facet_var = factor(facet_var, levels = facet_levels_B))

p_panelB <- ggplot(coverage_fig5, aes(x = TAD - 0.5, y = DRUG_KEY, fill = LESIONNAME2)) +
  geom_point(data = dplyr::filter(coverage_fig5, !is.na(TAD)),
             shape = 22, size = 1.8, stroke = 0.1) +
  facet_wrap(~facet_var, nrow = 2, ncol = 7, scales = "free_y", dir = "v",
             drop = FALSE) +
  scale_y_discrete(labels = function(x) sub("___.*", "", x),
                   expand = expansion(mult = c(0.1, 0.5))) +
  scale_fill_manual(values = c("Easy-to-treat" = blue, "Hard-to-treat" = orange),
                    guide = "none") +
  scale_x_continuous(breaks = seq(0, 24, 12)) +
  coord_cartesian(xlim = c(-0.5, 24)) +
  labs(x = "Predicted PK-PD drug coverage over one steady-state day (hours)",
       y = "Drugs\nin regimen") +
  theme_bw(base_size = 8) +
  theme(strip.text       = element_blank(),
        panel.grid       = element_blank(),
        axis.text        = element_text(size = 6),
        axis.title       = element_text(size = 8),
        panel.background = element_blank(),
        plot.margin      = margin(2, 5, 5, 5))

# ---------------------------------------------------------------------------
# Combine Panel A (KM) + Panel B (coverage) and save
# ---------------------------------------------------------------------------

p_figure5_combined <- plot_grid(
  p_figure5, p_panelB,
  ncol = 1, rel_heights = c(2.5, 1.5),
  labels = c("A.", "B."), label_size = 12,
  align = "v", axis = "lr"
)

fig5_path <- file.path(paths$results, "Fig5_clinical_outcome.png")
ggsave(fig5_path, p_figure5_combined, width = 14, height = 7.5, dpi = 300)
message("Saved: ", fig5_path)

sim_ppt_add_slide(x = p_figure5_combined, width = 14, height = 7.5)
message("Figure 5 appended to PPTX")
