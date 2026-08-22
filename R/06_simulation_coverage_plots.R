# Expects globals from 01_setup.R: lesion_level, drug_list_val, %notin%, paths, ppt_add_slide
# Expects globals from 02_import_data.R: drug_potency_MRT2_edited
# - coverage_plot_all_NO_METAB  (parent compounds)
# - coverage_plot_all_METAB     (metabolites only)

suppressPackageStartupMessages({
  library(xpose4)
})

sim_ppt_path <- file.path(paths$results, "simulation_coverage_figures.pptx")
.sim_ppt_append <- FALSE

sim_ppt_add_slide <- function(x = NULL, fun = NULL, width, height, ...) {
  if (!requireNamespace("export", quietly = TRUE)) {
    return(invisible(NULL))
  }
  if (is.null(x) && is.null(fun)) {
    stop("sim_ppt_add_slide: provide either x (ggplot) or fun (drawing function).")
  }
  if (!is.null(x) && !is.null(fun)) {
    stop("sim_ppt_add_slide: pass only one of x or fun.")
  }
  export::graph2ppt(
    x = x,
    fun = fun,
    file = sim_ppt_path,
    width = width,
    height = height,
    append = .sim_ppt_append,
    ...
  )
  .sim_ppt_append <<- TRUE
  invisible(sim_ppt_path)
}

# ---------------------------------------------------------------------------
# Load simulation data
# ---------------------------------------------------------------------------

sim_data_path <- here("data", "simulation_datasets", "20260802_simulated_data.RData")
if (!file.exists(sim_data_path)) {
  stop("Missing: ", sim_data_path, "\nSee data/simulation_datasets/")
}
load(sim_data_path)

clin_models <- here("data", "simulation_datasets", "clinical_pop_PK_models")

# ---------------------------------------------------------------------------
# GTX (PKPDsim CSV)
# ---------------------------------------------------------------------------

gtx_csv <- file.path(clin_models, "GTX", "df_gtx_sim_output.csv")
if (!file.exists(gtx_csv)) stop("Missing: ", gtx_csv)

df_gtx_conc <- read_csv(gtx_csv, show_col_types = FALSE) %>%
  rename(TIME = t) %>%
  mutate(CPLASMA = CP, ID = 1, DOSE = 400) %>%
  select(ID, DOSE, TIME, TAD, CP, CPLASMA, CLUNG, CCELL, CCAS) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "TAD", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("CPLASMA", CPname) ~ "Plasma",
    grepl("CLUNG", CPname)   ~ "Lung",
    grepl("CCELL", CPname)   ~ "Cellular lesion",
    grepl("CCAS", CPname)    ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "GTX",
         DRUG_name_long = "Gatifloxacin", FLAG = 1, II = "QD") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

# ---------------------------------------------------------------------------
# NONMEM simulations (CFZ; BDQ, TBAJ-587, TBAJ-876 + metabolites)
# ---------------------------------------------------------------------------

cfz_sdtab     <- file.path(clin_models, "CFZ",      "Final_Sim", "sdtab1")
bdq_sdtab     <- file.path(clin_models, "BDQ",      "Final_Sim", "sdtab663")
tbaj587_sdtab <- file.path(clin_models, "TBAJ-587", "Final_Sim", "sdtab7")
tbaj876_sdtab <- file.path(clin_models, "TBAJ-876", "Final_Sim", "sdtab0547")

for (f in c(cfz_sdtab, bdq_sdtab, tbaj587_sdtab, tbaj876_sdtab)) {
  if (!file.exists(f)) stop("Missing NONMEM table: ", f)
}

df_cfz_sim_output     <- read_nm_table(cfz_sdtab)
df_bdq_sim_output     <- read_nm_table(bdq_sdtab)
df_tbaj587_sim_output <- read_nm_table(tbaj587_sdtab)
df_tbaj876_sim_output <- read_nm_table(tbaj876_sdtab)

# --- Parent compounds ---

df_cfz_conc <- df_cfz_sim_output %>%
  filter(EVID == 0) %>%
  select(ID, DOSE, TIME, starts_with("C")) %>%
  select(-CWRES) %>%
  mutate(CPLASMA = CP) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("PLASMA", CPname) ~ "Plasma",
    grepl("LUNG",   CPname) ~ "Lung",
    grepl("CELL",   CPname) ~ "Cellular lesion",
    grepl("CASEUM", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "CFZ",
         DRUG_name_long = "Clofazimine", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

df_bdq_conc <- df_bdq_sim_output %>%
  filter(FLAG == 1, EVID == 0) %>%
  mutate(CP = CPplasma) %>%
  select(ID, DOSE, TIME, starts_with("CP")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("plasma", CPname) ~ "Plasma",
    grepl("lung",   CPname) ~ "Lung",
    grepl("cell",   CPname) ~ "Cellular lesion",
    grepl("caseum", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "BDQ",
         DRUG_name_long = "Bedaquiline", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

df_tbaj587_conc <- df_tbaj587_sim_output %>%
  filter(METAB == 1, EVID == 0) %>%
  mutate(DOSE = round(DOSE * 614.5, digits = 0)) %>%
  mutate(CP = CPplasma) %>%
  select(ID, DOSE, TIME, starts_with("CP")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("plasma", CPname) ~ "Plasma",
    grepl("lung",   CPname) ~ "Lung",
    grepl("cell",   CPname) ~ "Cellular lesion",
    grepl("caseum", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "TBAJ-587",
         DRUG_name_long = "TBAJ-587", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

df_tbaj876_conc <- df_tbaj876_sim_output %>%
  filter(METAB == 1, EVID == 0) %>%
  mutate(CP = Cparent) %>%
  select(ID, DOSE, TIME, CP, starts_with("Cp")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    CPname == "Cparent"      ~ "Plasma",
    grepl("lung",   CPname)  ~ "Lung",
    grepl("cell",   CPname)  ~ "Cellular lesion",
    grepl("caseum", CPname)  ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "TBAJ-876",
         DRUG_name_long = "TBAJ-876", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

# --- Metabolites ---

df_bdqm2_conc <- df_bdq_sim_output %>%
  filter(FLAG == 2, EVID == 0) %>%
  mutate(CP = CM2plasma) %>%
  select(ID, DOSE, TIME, CP, starts_with("CM2")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("plasma", CPname) ~ "Plasma",
    grepl("lung",   CPname) ~ "Lung",
    grepl("cell",   CPname) ~ "Cellular lesion",
    grepl("caseum", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "BDQ-M2",
         DRUG_name_long = "Bedaquiline-M2", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

df_tbaj587m3_conc <- df_tbaj587_sim_output %>%
  filter(METAB == 2, EVID == 0) %>%
  mutate(DOSE = round(DOSE * 614.5, digits = 0)) %>%
  mutate(CP = CM3plasma) %>%
  select(ID, DOSE, TIME, CP, starts_with("CM3")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    grepl("plasma", CPname) ~ "Plasma",
    grepl("lung",   CPname) ~ "Lung",
    grepl("cell",   CPname) ~ "Cellular lesion",
    grepl("caseum", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "TBAJ-587-M3",
         DRUG_name_long = "TBAJ-587-M3", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

df_tbaj876m3_conc <- df_tbaj876_sim_output %>%
  filter(METAB == 2, EVID == 0) %>%
  mutate(CP = Cmetab) %>%
  select(ID, DOSE, TIME, CP, starts_with("Cmeta")) %>%
  pivot_longer(cols = -c("ID", "DOSE", "TIME", "CP"),
               names_to = "CPname", values_to = "CLESION") %>%
  mutate(LESIONNAME = case_when(
    CPname == "Cmetab"      ~ "Plasma",
    grepl("lung",   CPname) ~ "Lung",
    grepl("cell",   CPname) ~ "Cellular lesion",
    grepl("caseum", CPname) ~ "Caseum"
  )) %>%
  mutate(SPECIES = "Corrected \nRabbit", DRUG = "TBAJ-876-M3",
         DRUG_name_long = "TBAJ-876-M3", FLAG = 1, II = "24") %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", II, sep = "_")) %>%
  filter(between(TIME, max(TIME) - 24, max(TIME))) %>%
  mutate(TAD = TIME - (max(TIME) - 24)) %>%
  select(ROW_ID, DRUG, DOSE, II, SPECIES, LESIONNAME, TIME, CP, CLESION, TAD, FLAG)

# ---------------------------------------------------------------------------
# Build sim_df: Replace PKPDsim CFZ/BDQ/TBAJ/GTX with NONMEM equivalents
# ---------------------------------------------------------------------------

sim_df <- simulation_collect_all_drugs %>%
  filter(!grepl("CFZ", DRUG)) %>%
  filter(!grepl("BDQ", DRUG)) %>%
  filter(!grepl("TBAJ", DRUG)) %>%
  filter(!grepl("GTX", DRUG)) %>%
  bind_rows(df_cfz_conc, 
            df_bdq_conc, df_tbaj587_conc, df_tbaj876_conc,
            df_bdqm2_conc, df_tbaj587m3_conc, df_tbaj876m3_conc,
            df_gtx_conc) %>%
  mutate(DRUG = ifelse(DRUG == "INH_s" | DRUG == "INH_f", "INH", DRUG)) %>%
  mutate(FLAG = ifelse(DRUG %in% drug_list_val & SPECIES == "Human", 1,
                       ifelse(DRUG %notin% drug_list_val & grepl("Corrected", SPECIES), 1, 0))) %>%
  filter(FLAG == 1) %>%
  mutate(SPECIES = ifelse(SPECIES == "Corrected \nRabbit", "Corrected Rabbit", SPECIES)) %>%
  filter(LESIONNAME %in% c("Plasma", "Lung", "Cellular lesion", "Caseum")) %>%
  mutate(DRUG = case_when(
    DRUG == "OPC-167832" ~ "QBS",
    DRUG == "QBD"        ~ "QBS",
    DRUG == "GSK-656"    ~ "GFB",
    DRUG == "mCLB-073"   ~ "TBD-11",
    .default = DRUG
  ))

target_labels <- drug_potency_MRT2_edited %>%
  select(DRUG, DRUG_name_long, LESIONNAME, contains("mgL")) %>%
  mutate(DRUG = case_when(
    DRUG == "OPC-167832" ~ "QBS",
    DRUG == "GSK-656"    ~ "GFB",
    DRUG == "mCLB-073"   ~ "TBD-11",
    .default = DRUG
  )) %>%
  pivot_longer(cols = contains("mgL"), names_to = "LABEL", values_to = "TARGET")

sim_df <- sim_df %>%
  left_join(target_labels, relationship = "many-to-many") %>%
  mutate(LESIONNAME = factor(LESIONNAME, levels = lesion_level))

message("sim_df drugs: ", paste(sort(unique(sim_df$DRUG)), collapse = ", "))

# ---------------------------------------------------------------------------
# Coverage computation
# ---------------------------------------------------------------------------

drug_classes <- list(
  "AMK" = "aminoglycoside", "BTZ-043" = "zother", "CFZ" = "zother",
  "DLM" = "nitroimidazole", "EMB" = "zother", "GSK-656" = "zother",
  "GFB" = "zother", "INH" = "zother", "KAN" = "aminoglycoside",
  "LZD" = "oxazolidinone", "mCLB-073" = "zother", "TBD-11" = "zother",
  "MXF" = "fluoroquinolone", "OPC-167832" = "zother", "QBS" = "zother",
  "PMD" = "nitroimidazole", "PZA" = "zother", "RBT" = "rifamycin",
  "RIF" = "rifamycin", "RPT" = "rifamycin",
  "SZD" = "oxazolidinone", "SZD-M1" = "oxazolidinone", "TBI-223" = "oxazolidinone",
  "BDQ" = "diarylquinoline", "TBAJ-587" = "diarylquinoline", "TBAJ-876" = "diarylquinoline",
  "BDQ-M2" = "diarylquinoline", "TBAJ-587-M3" = "diarylquinoline",
  "TBAJ-876-M3" = "diarylquinoline", "GTX" = "fluoroquinolone"
)

coverage_all_drugs_intermediate <- sim_df %>%
  mutate(II_num = ifelse(II == "BID", 12, 24)) %>%
  mutate(ROW_ID = paste(DRUG, DOSE, "mg", sep = " ")) %>%
  mutate(LABEL = ifelse(LABEL == "casMBC90_mgL", "FIG4_casMBC90_mgL", LABEL)) %>%
  filter(LABEL %in% c("MRT_mgL", "FIG4_casMBC90_mgL")) %>%
  mutate(
    coverage_above_MRT      = ifelse(CLESION > TARGET & LABEL == "MRT_mgL", 1, 0),
    coverage_above_CasMBC90 = ifelse(CLESION > TARGET & LABEL == "FIG4_casMBC90_mgL", 1, 0)
  ) %>%
  mutate(CLASS = as.character(drug_classes[DRUG]))

coverage_all_drugs <- coverage_all_drugs_intermediate %>%
  filter(coverage_above_MRT == 1 | coverage_above_CasMBC90 == 1) %>%
  mutate(LESIONNAME = factor(LESIONNAME, levels = lesion_level)) %>%
  group_by(DRUG, DOSE) %>%
  mutate(ROW_ID = paste(CLASS, DRUG, DOSE, "mg", sep = " "))

# ---------------------------------------------------------------------------
# Shared ggplot theme for coverage plots
# ---------------------------------------------------------------------------

metab_drugs <- c("TBAJ-587-M3", "TBAJ-876-M3", "BDQ-M2", "SZD-M1")

coverage_theme <- theme_classic(base_size = 10) +
  theme(
    legend.background  = element_rect("transparent"),
    strip.background   = element_rect("transparent"),
    legend.text        = element_text(size = 5),
    legend.margin      = margin(0, 0, 0, 0.1),
    legend.direction   = "vertical",
    legend.key.height  = unit(0.005, "in"),
    legend.key.size    = unit(0.1, "in"),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid         = element_blank(),
    axis.title         = element_text(size = 7),
    axis.text          = element_text(size = 7),
    strip.text         = element_text(size = 8, margin = margin(0.1, 0, 0.05, 0, "cm")),
    plot.margin        = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
  )

coverage_fill <- scale_fill_manual(values = c(
  "white"       = "steelblue4",
  "steelblue4"  = "steelblue4",
  "slategray3"  = "slategray3"
))

lesion_header_levels <- c(
  "Plasma\nconcentrations over\nMIC",
  "Lung\nconcentrations over\nMacIC90",
  "Cellular lesion\nconcentrations over\nMacIC90",
  "Caseum\nconcentrations over\ncasMBC50,casMBC90"
)

add_lesion_header <- function(df) {
  df %>%
    mutate(LESIONNAME_header = case_when(
      LESIONNAME == "Plasma"          ~ "Plasma\nconcentrations over\nMIC",
      LESIONNAME == "Lung"            ~ "Lung\nconcentrations over\nMacIC90",
      LESIONNAME == "Cellular lesion" ~ "Cellular lesion\nconcentrations over\nMacIC90",
      LESIONNAME == "Caseum"          ~ "Caseum\nconcentrations over\ncasMBC50,casMBC90"
    )) %>%
    mutate(LESIONNAME_header = factor(LESIONNAME_header, levels = lesion_header_levels))
}

add_fill_color <- function(df) {
  df %>%
    mutate(fill_color = case_when(
      LESIONNAME_header == "Caseum\nconcentrations over\ncasMBC50,casMBC90" &
        coverage_above_CasMBC90 == 1 ~ "steelblue4",
      LESIONNAME_header == "Caseum\nconcentrations over\ncasMBC50,casMBC90" &
        coverage_above_MRT == 1 & coverage_above_CasMBC90 != 1 ~ "slategray3",
      coverage_above_MRT == 1 ~ "steelblue4",
      TRUE ~ "white"
    ))
}

# ---------------------------------------------------------------------------
# Helper: build one coverage plot
# ---------------------------------------------------------------------------

make_coverage_plot <- function(keep_metab) {
  if (keep_metab) {
    drug_filter <- function(df) filter(df, DRUG %in% metab_drugs)
  } else {
    drug_filter <- function(df) filter(df, DRUG %notin% metab_drugs)
  }

  blanks <- cross_join(
    coverage_all_drugs_intermediate %>%
      select(ROW_ID, CLASS, DRUG, DOSE) %>%
      unique(),
    expand_grid(
      TAD = seq(0, 24),
      LESIONNAME = factor(c("Plasma", "Lung", "Cellular lesion", "Caseum"),
                          levels = lesion_level)
    )
  ) %>%
    mutate(DRUG = case_when(
      DRUG == "OPC-167832" ~ "QBS",
      DRUG == "GSK-656"    ~ "GFB",
      DRUG == "mCLB-073"   ~ "TBD-11",
      .default = DRUG
    )) %>%
    drug_filter() %>%
    ungroup() %>%
    arrange(desc(CLASS), desc(DRUG), desc(DOSE)) %>%
    mutate(ROW_ID = paste(CLASS, DRUG, DOSE, "mg", sep = " "))

  coverage_all_drugs %>%
    drug_filter() %>%
    filter(LESIONNAME %notin% c("Cavity wall", "Caseous lesion")) %>%
    mutate(DRUG = case_when(
      DRUG == "OPC-167832" ~ "QBS",
      DRUG == "GSK-656"    ~ "GFB",
      DRUG == "mCLB-073"   ~ "TBD-11",
      TRUE ~ DRUG
    )) %>%
    mutate(DRUG = ifelse(DRUG == "INH_s", "INH", DRUG)) %>%
    add_lesion_header() %>%
    mutate(ROW_ID = paste(CLASS, DRUG, DOSE, "mg", sep = " ")) %>%
    ungroup() %>%
    arrange(desc(CLASS), desc(DRUG), desc(DOSE)) %>%
    mutate(ROW_ID = factor(ROW_ID, unique(ROW_ID))) %>%
    add_fill_color() %>%
    arrange(fill_color) %>%
    ggplot(aes(x = TAD + 0.5, y = ROW_ID, fill = fill_color)) +
    geom_point(
      data = blanks %>%
        ungroup() %>%
        arrange(desc(CLASS), desc(DRUG), desc(DOSE)) %>%
        mutate(ROW_ID = factor(ROW_ID, unique(ROW_ID))),
      aes(x = TAD + 0.5, y = ROW_ID),
      color = "black", shape = 22, size = 1.8, stroke = 0.2, inherit.aes = FALSE
    ) +
    geom_point(shape = 22, size = 1.8, stroke = 0.2) +
    coverage_theme +
    scale_x_continuous(breaks = seq(0, 24, 4), labels = seq(0, 24, 4), limits = c(0, 24)) +
    facet_wrap(~LESIONNAME_header, nrow = 1) +
    labs(y = "", x = "Time after dose (hours)") +
    coverage_fill
}

# ---------------------------------------------------------------------------
# Build & save plots
# ---------------------------------------------------------------------------

p_coverage_no_metab <- make_coverage_plot(keep_metab = FALSE)
p_coverage_metab    <- make_coverage_plot(keep_metab = TRUE)

ggsave(plot = p_coverage_no_metab,
       filename = file.path(paths$results, "coverage_plot_all_NO_METAB.png"),
       width = 7.4, height = 5, units = "in", dpi = 300)

ggsave(plot = p_coverage_metab,
       filename = file.path(paths$results, "coverage_plot_all_METAB.png"),
       width = 7.4, height = 2, units = "in", dpi = 300)

message("Saved: ", file.path(paths$results, "coverage_plot_all_NO_METAB.png"))
message("Saved: ", file.path(paths$results, "coverage_plot_all_METAB.png"))

sim_ppt_add_slide(x = p_coverage_no_metab, width = 7.4, height = 5)
sim_ppt_add_slide(x = p_coverage_metab,    width = 7.4, height = 2)

# ---------------------------------------------------------------------------
# Figure 4: Publication-quality coverage plot (NO METAB)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggtext)
  library(cowplot)
})

fig4_class_order <- c("aminoglycoside", "diarylquinoline", "fluoroquinolone",
                      "nitroimidazole", "oxazolidinone", "rifamycin", "zother")

fig4_class_labels <- c(
  aminoglycoside  = "aminoglycosides",
  diarylquinoline = "diarylquinolines",
  fluoroquinolone = "fluoroquinolones",
  nitroimidazole  = "nitroimidazoles",
  oxazolidinone   = "oxazolidinones",
  rifamycin       = "rifamycins",
  zother          = "other structural\nclasses"
)

fig4_class_colors <- c(
  aminoglycoside  = "#00BFC4",
  diarylquinoline = "#2B6CB0",
  fluoroquinolone = "#38A1DB",
  nitroimidazole  = "#E76BF3",
  oxazolidinone   = "#00BA38",
  rifamycin       = "#F8766D",
  zother          = "#7F7F7F"
)

fig4_headers <- c(
  "Plasma\nconcentrations over\nMIC",
  "Lung\nconcentrations over\nMacIC90",
  "Cellular lesion\nconcentrations over\nMacIC90",
  "Caseum\nconcentrations over\ncasMBC50, casMBC90"
)

# Row ordering: class -> drug (alpha) -> dose (ascending)
fig4_rows <- coverage_all_drugs_intermediate %>%
  filter(DRUG %notin% metab_drugs) %>%
  select(DRUG, DOSE, CLASS) %>%
  distinct() %>%
  mutate(DRUG_DISPLAY = case_when(DRUG == "QBS" ~ "QBD", TRUE ~ DRUG)) %>%
  mutate(class_rank = match(CLASS, fig4_class_order)) %>%
  arrange(class_rank, DRUG_DISPLAY, DOSE) %>%
  mutate(DISPLAY = paste(DRUG_DISPLAY, DOSE, "mg"))

fig4_levels <- rev(fig4_rows$DISPLAY)
fig4_rows$y_pos <- match(fig4_rows$DISPLAY, fig4_levels)

# Class bracket positions (numeric y coords matching discrete factor positions)
fig4_bracket_df <- fig4_rows %>%
  group_by(CLASS) %>%
  summarise(ymin = min(y_pos) - 0.4, ymax = max(y_pos) + 0.4,
            ymid = (min(y_pos) + max(y_pos)) / 2, .groups = "drop") %>%
  mutate(label = fig4_class_labels[CLASS],
         color = fig4_class_colors[CLASS],
         LESIONNAME_header = factor(fig4_headers[1], levels = fig4_headers))

# Colored y-axis labels via ggtext
fig4_label_colors <- setNames(fig4_class_colors[fig4_rows$CLASS], fig4_rows$DISPLAY)
fig4_colored_labels <- setNames(
  sprintf("<span style='color:%s'>%s</span>",
          fig4_label_colors[fig4_levels], fig4_levels),
  fig4_levels
)

# Background grid
fig4_blanks <- cross_join(
  fig4_rows %>% select(DRUG_DISPLAY, DOSE),
  expand_grid(TAD = seq(0, 24),
              LESIONNAME = factor(c("Plasma", "Lung", "Cellular lesion", "Caseum"),
                                  levels = lesion_level))
) %>%
  mutate(DISPLAY = factor(paste(DRUG_DISPLAY, DOSE, "mg"), levels = fig4_levels))

# Coverage data
fig4_cov <- coverage_all_drugs %>%
  ungroup() %>%
  filter(DRUG %notin% metab_drugs) %>%
  filter(LESIONNAME %notin% c("Cavity wall", "Caseous lesion")) %>%
  # casMBC90 is a caseum-specific target; keep those rows only for the Caseum panel
  # so Plasma/Lung/Cellular render strictly over MIC/macIC90 (their MRT_mgL target)
  filter(!(LABEL == "FIG4_casMBC90_mgL" & LESIONNAME != "Caseum")) %>%
  mutate(DRUG_DISPLAY = case_when(DRUG == "QBS" ~ "QBD", TRUE ~ DRUG),
         DRUG = ifelse(DRUG == "INH_s", "INH", DRUG)) %>%
  mutate(LESIONNAME_header = factor(case_when(
    LESIONNAME == "Plasma"          ~ fig4_headers[1],
    LESIONNAME == "Lung"            ~ fig4_headers[2],
    LESIONNAME == "Cellular lesion" ~ fig4_headers[3],
    LESIONNAME == "Caseum"          ~ fig4_headers[4]
  ), levels = fig4_headers)) %>%
  mutate(fill_color = case_when(
    LESIONNAME_header == fig4_headers[4] & coverage_above_CasMBC90 == 1 ~ "steelblue4",
    LESIONNAME_header == fig4_headers[4] & coverage_above_MRT == 1 &
      coverage_above_CasMBC90 != 1 ~ "slategray3",
    coverage_above_MRT == 1 ~ "steelblue4",
    TRUE ~ "white"
  )) %>%
  mutate(DISPLAY = factor(paste(DRUG_DISPLAY, DOSE, "mg"), levels = fig4_levels)) %>%
  arrange(fill_color)

# Build Figure 4 using patchwork: bracket panel (left) + coverage (right)

# --- Main coverage panel ---
p_fig4_main <- ggplot(fig4_cov, aes(x = TAD + 0.5, y = DISPLAY, fill = fill_color)) +
  geom_point(data = fig4_blanks, aes(x = TAD + 0.5, y = DISPLAY),
             color = "black", shape = 22, size = 1.8, stroke = 0.2,
             inherit.aes = FALSE) +
  geom_point(shape = 22, size = 1.8, stroke = 0.2) +
  scale_fill_manual(values = c("white" = "steelblue4", "steelblue4" = "steelblue4",
                                "slategray3" = "slategray3")) +
  scale_x_continuous(breaks = seq(0, 24, 4), labels = seq(0, 24, 4), limits = c(0, 24)) +
  scale_y_discrete(labels = fig4_colored_labels) +
  facet_wrap(~LESIONNAME_header, nrow = 1) +
  labs(y = "", x = "Predicted PK-PD drug coverage over one steady-state day (hours)") +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "none",
    panel.grid         = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect("transparent"),
    axis.text.y        = element_markdown(size = 6),
    axis.text.x        = element_text(size = 7),
    axis.title         = element_text(size = 7),
    strip.text         = element_text(size = 8, margin = margin(0.1, 0, 0.05, 0, "cm")),
    plot.margin        = margin(5, 5, 5, 2, "pt")
  )

# --- Class bracket panel (same discrete y-axis for alignment) ---
fig4_bracket_plot_df <- fig4_bracket_df %>% select(-LESIONNAME_header)

p_fig4_brackets <- ggplot(data.frame(DISPLAY = factor(fig4_levels, levels = fig4_levels)),
                           aes(y = DISPLAY)) +
  geom_blank(aes(x = 1)) +
  geom_segment(data = fig4_bracket_plot_df,
               aes(x = 0.85, xend = 0.85, y = ymin, yend = ymax),
               color = fig4_bracket_plot_df$color, linewidth = 1.5,
               inherit.aes = FALSE) +
  geom_text(data = fig4_bracket_plot_df,
            aes(x = 0.75, y = ymid, label = label),
            color = fig4_bracket_plot_df$color,
            size = 2.3, hjust = 1, fontface = "italic", lineheight = 0.85,
            inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-0.2, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    axis.text.y  = element_blank(),
    plot.margin  = margin(5, 0, 5, 2, "pt")
  )

# --- Combine via cowplot (horizontal alignment keeps y-axes matching) ---
p_figure4 <- cowplot::plot_grid(p_fig4_brackets, p_fig4_main,
                                nrow = 1, rel_widths = c(1.2, 8),
                                align = "h", axis = "tb")

ggsave(plot = p_figure4,
       filename = file.path(paths$results, "Fig4_coverage_NO_METAB.png"),
       width = 8.5, height = 5, units = "in", dpi = 300)

message("Saved: ", file.path(paths$results, "Fig4_coverage_NO_METAB.png"))

sim_ppt_add_slide(x = p_figure4, width = 8.5, height = 5)

# ---------------------------------------------------------------------------
# Figure S9: Publication-quality coverage plot (METAB only)
# ---------------------------------------------------------------------------

figS9_rows <- coverage_all_drugs_intermediate %>%
  filter(DRUG %in% metab_drugs) %>%
  select(DRUG, DOSE, CLASS) %>%
  distinct() %>%
  mutate(DRUG_DISPLAY = DRUG) %>%
  mutate(class_rank = match(CLASS, fig4_class_order)) %>%
  arrange(class_rank, DRUG_DISPLAY, DOSE) %>%
  mutate(DISPLAY = paste(DRUG_DISPLAY, DOSE, "mg"))

figS9_levels <- rev(figS9_rows$DISPLAY)
figS9_rows$y_pos <- match(figS9_rows$DISPLAY, figS9_levels)

figS9_bracket_df <- figS9_rows %>%
  group_by(CLASS) %>%
  summarise(ymin = min(y_pos) - 0.4, ymax = max(y_pos) + 0.4,
            ymid = (min(y_pos) + max(y_pos)) / 2, .groups = "drop") %>%
  mutate(label = fig4_class_labels[CLASS],
         color = fig4_class_colors[CLASS])

figS9_label_colors <- setNames(fig4_class_colors[figS9_rows$CLASS], figS9_rows$DISPLAY)
figS9_colored_labels <- setNames(
  sprintf("<span style='color:%s'>%s</span>",
          figS9_label_colors[figS9_levels], figS9_levels),
  figS9_levels
)

figS9_blanks <- cross_join(
  figS9_rows %>% select(DRUG_DISPLAY, DOSE),
  expand_grid(TAD = seq(0, 24),
              LESIONNAME = factor(c("Plasma", "Lung", "Cellular lesion", "Caseum"),
                                  levels = lesion_level))
) %>%
  mutate(DISPLAY = factor(paste(DRUG_DISPLAY, DOSE, "mg"), levels = figS9_levels))

figS9_cov <- coverage_all_drugs %>%
  ungroup() %>%
  filter(DRUG %in% metab_drugs) %>%
  filter(LESIONNAME %notin% c("Cavity wall", "Caseous lesion")) %>%
  # casMBC90 is a caseum-specific target; keep those rows only for the Caseum panel
  filter(!(LABEL == "FIG4_casMBC90_mgL" & LESIONNAME != "Caseum")) %>%
  mutate(DRUG_DISPLAY = DRUG) %>%
  mutate(LESIONNAME_header = factor(case_when(
    LESIONNAME == "Plasma"          ~ fig4_headers[1],
    LESIONNAME == "Lung"            ~ fig4_headers[2],
    LESIONNAME == "Cellular lesion" ~ fig4_headers[3],
    LESIONNAME == "Caseum"          ~ fig4_headers[4]
  ), levels = fig4_headers)) %>%
  mutate(fill_color = case_when(
    LESIONNAME_header == fig4_headers[4] & coverage_above_CasMBC90 == 1 ~ "steelblue4",
    LESIONNAME_header == fig4_headers[4] & coverage_above_MRT == 1 &
      coverage_above_CasMBC90 != 1 ~ "slategray3",
    coverage_above_MRT == 1 ~ "steelblue4",
    TRUE ~ "white"
  )) %>%
  mutate(DISPLAY = factor(paste(DRUG_DISPLAY, DOSE, "mg"), levels = figS9_levels)) %>%
  arrange(fill_color)

# --- Main coverage panel ---
p_figS9_main <- ggplot(figS9_cov, aes(x = TAD + 0.5, y = DISPLAY, fill = fill_color)) +
  geom_point(data = figS9_blanks, aes(x = TAD + 0.5, y = DISPLAY),
             color = "black", shape = 22, size = 1.8, stroke = 0.2,
             inherit.aes = FALSE) +
  geom_point(shape = 22, size = 1.8, stroke = 0.2) +
  scale_fill_manual(values = c("white" = "steelblue4", "steelblue4" = "steelblue4",
                                "slategray3" = "slategray3")) +
  scale_x_continuous(breaks = seq(0, 24, 4), labels = seq(0, 24, 4), limits = c(0, 24)) +
  scale_y_discrete(labels = figS9_colored_labels) +
  facet_wrap(~LESIONNAME_header, nrow = 1) +
  labs(y = "", x = "Predicted PK-PD drug coverage over one steady-state day (hours)") +
  theme_classic(base_size = 10) +
  theme(
    legend.position    = "none",
    panel.grid         = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect("transparent"),
    axis.text.y        = element_markdown(size = 6),
    axis.text.x        = element_text(size = 7),
    axis.title         = element_text(size = 7),
    strip.text         = element_text(size = 8, margin = margin(0.1, 0, 0.05, 0, "cm")),
    plot.margin        = margin(5, 5, 5, 2, "pt")
  )

# --- Class bracket panel ---
p_figS9_brackets <- ggplot(data.frame(DISPLAY = factor(figS9_levels, levels = figS9_levels)),
                            aes(y = DISPLAY)) +
  geom_blank(aes(x = 1)) +
  geom_segment(data = figS9_bracket_df,
               aes(x = 0.85, xend = 0.85, y = ymin, yend = ymax),
               color = figS9_bracket_df$color, linewidth = 1.5,
               inherit.aes = FALSE) +
  geom_text(data = figS9_bracket_df,
            aes(x = 0.75, y = ymid, label = label),
            color = figS9_bracket_df$color,
            size = 2.3, hjust = 1, fontface = "italic", lineheight = 0.85,
            inherit.aes = FALSE) +
  scale_x_continuous(limits = c(-0.2, 1), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    axis.text.y  = element_blank(),
    plot.margin  = margin(5, 0, 5, 2, "pt")
  )

# --- Combine ---
p_figureS9 <- cowplot::plot_grid(p_figS9_brackets, p_figS9_main,
                                  nrow = 1, rel_widths = c(1.2, 8),
                                  align = "h", axis = "tb")

ggsave(plot = p_figureS9,
       filename = file.path(paths$results, "FigS9_coverage_METAB.png"),
       width = 8.5, height = 2.5, units = "in", dpi = 300)

message("Saved: ", file.path(paths$results, "FigS9_coverage_METAB.png"))

sim_ppt_add_slide(x = p_figureS9, width = 8.5, height = 2.5)
