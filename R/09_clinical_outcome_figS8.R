# Self-contained: only requires here() paths and 3 CSV data files
# Expects sim_ppt_add_slide() to be in the global environment (from 06)

# Same 7 facets as Figure 5
# Stratified by cavitary status only (contrast with Figure 5 composite stratification)

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
# Risk stratification — cavitation only
# ---------------------------------------------------------------------------

DF_TBR <- DF_TBR %>%
  filter(CAVITY != -99, SMEAR != -99) %>%
  mutate(GROUP = ifelse(CAVITY == 1, "Hard", "Moderate/easy"),
         LUNGGRADING1 = -99)

DF_S31 <- DF_S31 %>%
  filter(XPERT != -99, CXREXTNT != -99) %>%
  mutate(GROUP = ifelse(CAVITY == 1, "Hard", "Moderate/easy"),
         LUNGGRADING1 = -99)

DF_RIFASHORT <- DF_RIFASHORT %>%
  mutate(LUNGGRADING  = as.factor(LUNGGRADING),
         LUNGGRADING1 = as.factor(ifelse(is.na(LUNGGRADING), 1,
                                         ifelse(LUNGGRADING == 4, 1, 0))))

DF_RIFASHORT_XPERT <- DF_RIFASHORT %>%
  filter(XPERT != -99) %>%
  mutate(GROUP = ifelse(CAVITY == 1, "Hard", "Moderate/easy"))

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
# BPaL (Nix/ZeNix) — cavitation-only stratification
# ---------------------------------------------------------------------------

all_mdr <- read.csv(file.path(clin_data, "All_validationMDR.csv"), header = TRUE)
all_mdr$DV <- as.numeric(as.character(all_mdr$DV))
all_mdr <- all_mdr %>%
  rename(TIMEM = TIME2) %>%
  mutate(RISK_MDR3 = ifelse(CAVITY < 2, "Low", "Moderate/High"))

DF_BPAL <- all_mdr %>%
  group_by(SUBJID) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(TIME = TIMEM) %>%
  rename(Status = DV) %>%
  mutate(RISK_GROUP = ifelse(RISK_MDR3 == "Moderate/High", "Hard", "Moderate/easy")) %>%
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
  legend.labs   = c("Minimal cavitary disease", "Extensive cavitary disease"),
  break.time.by = 6
)

ALL_plot <- ALL_plot +
  coord_cartesian(xlim = c(0, 18), ylim = c(0.75, 1)) +
  labs(x = "Months since treatment initiation",
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

p_figureS8 <- ALL_plot +
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
# Save
# ---------------------------------------------------------------------------

figS8_path <- file.path(paths$results, "FigS8_clinical_outcome_cavitation.png")
ggsave(figS8_path, p_figureS8, width = 14, height = 4.5, dpi = 300)
message("Saved: ", figS8_path)

sim_ppt_add_slide(x = p_figureS8, width = 14, height = 4.5)
message("Figure S8 appended to PPTX")
