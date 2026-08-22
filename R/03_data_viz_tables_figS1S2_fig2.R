table1_rabbit <- left_join(
  all_rabbit_obs %>%
    distinct(DRUG, CID) %>%
    group_by(DRUG) %>%
    summarise(Subject_N = n(), .groups = "drop"),
  all_rabbit_obs %>%
    filter(LESIONGROUP == "Plasma") %>%
    count(DRUG, LESIONGROUP) %>%
    pivot_wider(names_from = LESIONGROUP, values_from = n)
) %>%
  left_join(
    all_rabbit_obs %>%
      filter(LESIONGROUP != "Plasma") %>%
      mutate(LESIONGROUP = "Lesion") %>%
      count(DRUG, LESIONGROUP) %>%
      pivot_wider(names_from = LESIONGROUP, values_from = n),
    by = "DRUG"
  )

table1_human <- left_join(
  all_human_obs %>%
    distinct(DRUG, CID) %>%
    group_by(DRUG) %>%
    summarise(Subject_N = n(), .groups = "drop"),
  all_human_obs %>%
    filter(LESIONGROUP == "Plasma") %>%
    count(DRUG, LESIONGROUP) %>%
    pivot_wider(names_from = LESIONGROUP, values_from = n)
) %>%
  left_join(
    all_human_obs %>%
      filter(LESIONGROUP != "Plasma") %>%
      mutate(LESIONGROUP = "Lesion") %>%
      count(DRUG, LESIONGROUP) %>%
      pivot_wider(names_from = LESIONGROUP, values_from = n),
    by = "DRUG"
  )

table1_combine <- bind_rows(
  table1_human %>% mutate(Species = "Human", .before = DRUG),
  table1_rabbit %>% mutate(Species = "Rabbit", .before = DRUG)
)

write_csv(table1_combine, file.path(paths$results, "TableS1_data_available.csv"))

drug_potency_values <- drug_potency_MRT2_edited %>%
  select(-any_of(c("LESIONNAME", "MRT_mgL"))) %>%
  distinct() %>%
  as_tibble()

write_csv(drug_potency_values, file.path(paths$results, "TableS4_in_vitro_data.csv"))

all_raw_data_dict <- tribble(
  ~Name, ~Description,
  "DRUG", "Drug name",
  "CID", "Individual identificaiton number",
  "TAD", "time after dose in hours",
  "DV", "Concentration in mg/L (dependent variable)",
  "LESIONNAME", "The original label of the sample provided by HMH/CDI",
  "LESIONGROUP", "A grouping method",
  "METHOD", "Method for drug quantification; LCMS or LCM",
  "SPECIES", "Species; human or rabbit"
)
write_csv(all_raw_data_dict, file.path(paths$source_data, "all_raw_data_pk_data_dictionary.csv"))

Human_PK_plot_1A <- lapply(sort(unique(all_human_obs$DRUG)), function(drug) {
  all_human_obs %>%
    filter(DRUG == drug, DV > 0) %>%
    mutate(LESIONGROUP = ifelse(LESIONGROUP == "Cellular", "Cellular Lesion", LESIONGROUP)) %>%
    filter(LESIONGROUP %in% c("Plasma", "Lung", "Cellular Lesion", "Caseum")) %>%
    mutate(LESIONGROUP = factor(LESIONGROUP, levels = c("Plasma", "Lung", "Cellular Lesion", "Caseum"))) %>%
    ggplot(mapping = aes(x = TIME, y = DV, color = METHOD, group = METHOD)) +
    geom_point(size = 1, alpha = 0.2) +
    stat_summary_bin(geom = "line", fun = "median", breaks = c(-0.1, 4, 8, 12, 20, 24, 38), size = 0.8) +
    facet_wrap(~LESIONGROUP, nrow = 1) +
    ggtitle(drug) +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(0, 36, 6), limits = c(0, NA)) +
    scale_color_manual(values = c(COLOR_map[2], COLOR_map[1])) +
    theme_bw(base_size = 10) +
    theme(
      legend.background = element_rect("transparent"),
      strip.background = element_rect("transparent"),
      legend.text = element_text(size = 7),
      legend.margin = margin(0, 0, 0, 0),
      legend.direction = "vertical",
      legend.key.height = unit(0.005, "in"),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      axis.title = element_blank(),
      strip.text = element_text(size = 8, margin = margin(0.1, 0, 0.05, 0, "cm")),
      plot.title = element_text(size = 9, margin = margin(0.1, 0, 0.05, 0, "cm")),
      plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
    )
})

tiff(file.path(paths$results, "FigS2.tiff"), width = 6.5, height = 8, units = "in", res = 300, compression = "lzw")
grid.arrange(
  ggarrange(plotlist = Human_PK_plot_1A, nrow = 7, common.legend = TRUE, legend = "right"),
  bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 10)),
  left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 10), rot = 90)
)
dev.off()
ppt_add_slide(
  fun = function() {
    grid.arrange(
      ggarrange(plotlist = Human_PK_plot_1A, nrow = 7, common.legend = TRUE, legend = "right"),
      bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 10)),
      left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 10), rot = 90)
    )
  },
  width = 6.5,
  height = 8
)

both_species_slim_df_raw_data <- bind_rows(
  all_human_obs %>%
    select(DRUG, CID, TAD, DV, LESIONNAME, LESIONGROUP, METHOD) %>%
    mutate(SPECIES = "Clinical"),
  all_rabbit_obs %>%
    filter(SET == "Validation") %>%
    select(DRUG, CID, TAD, DV, LESIONNAME, LESIONGROUP, METHOD) %>%
    mutate(METHOD = ifelse(METHOD == 0, "LCMS", "LCM"), SPECIES = "Rabbit")
) %>%
  mutate(METHOD = ifelse(LESIONGROUP == "Plasma", "LCMS", METHOD))

Both_PK_plot_1 <- lapply(sort(unique(both_species_slim_df_raw_data$DRUG)), function(drug) {
  both_species_slim_df_raw_data %>%
    filter(DRUG == drug, DV > 0) %>%
    mutate(LESIONGROUP = ifelse(LESIONGROUP == "Cellular", "Cellular Lesion", LESIONGROUP)) %>%
    filter(LESIONGROUP %in% c("Plasma", "Lung", "Cellular Lesion", "Caseum")) %>%
    mutate(LESIONGROUP = factor(LESIONGROUP, levels = c("Plasma", "Lung", "Cellular Lesion", "Caseum"))) %>%
    ggplot(mapping = aes(x = TAD, y = DV, color = SPECIES, group = SPECIES)) +
    geom_point(size = 1, alpha = 0.2) +
    stat_summary_bin(geom = "line", fun = "median", breaks = c(-0.1, 4, 8, 12, 20, 24, 38, 54, 90), size = 0.8) +
    facet_wrap(~LESIONGROUP, nrow = 1) +
    ggtitle(drug) +
    scale_y_log10() +
    scale_x_continuous(
      breaks = if (identical(drug, "GFB")) c(0, 24, 48) else seq(0, 72, 12),
      limits = c(0, NA)
    ) +
    scale_color_manual(values = c("#FFA500", "#0072C6")) +
    theme_bw(base_size = 9) +
    theme(
      legend.background = element_rect("transparent"),
      strip.background = element_rect("transparent"),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.1, "in"),
      legend.margin = margin(0, 0, 0, 0.01),
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = 9),
      axis.title = element_blank(),
      strip.text = element_text(size = 8, margin = margin(0.1, 0, 0.05, 0, "cm")),
      plot.title = element_text(size = 9, margin = margin(0.1, 0, 0.05, 0, "cm")),
      plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
    )
})

tiff(file.path(paths$results, "FigS1.tiff"), width = 6.5, height = 8, units = "in", res = 300, compression = "lzw")
grid.arrange(
  ggarrange(plotlist = Both_PK_plot_1, nrow = 7, common.legend = TRUE, legend = "right"),
  bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 10)),
  left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 10), rot = 90)
)
dev.off()
ppt_add_slide(
  fun = function() {
    grid.arrange(
      ggarrange(plotlist = Both_PK_plot_1, nrow = 7, common.legend = TRUE, legend = "right"),
      bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 10)),
      left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 10), rot = 90)
    )
  },
  width = 6.5,
  height = 8
)

all_raw_data <- bind_rows(
  all_human_obs %>%
    select(DRUG, CID, TAD, DV, LESIONNAME, LESIONGROUP, METHOD) %>%
    mutate(SPECIES = "XClinical"),
  all_rabbit_obs %>%
    select(DRUG, CID, TAD, DV, LESIONNAME, LESIONGROUP, METHOD) %>%
    mutate(METHOD = ifelse(METHOD == 0, "LCMS", "LCM"), SPECIES = "Rabbit")
) %>%
  mutate(
    DRUG = ifelse(DRUG == "OPC-167832", "QBS", DRUG),
    DRUG = ifelse(DRUG == "mCLB-073", "TBD-11", DRUG),
    DRUG = ifelse(DRUG == "GSK-656", "GFB", DRUG)
  ) %>%
  arrange(SPECIES) %>%
  mutate(METHOD = ifelse(LESIONGROUP == "Plasma", "LCMS", METHOD))

write_csv(all_raw_data, file.path(paths$source_data, "all_raw_data_pk.csv"))

SPECIES_color_map <- c("XClinical" = "#FFA500", "Rabbit" = "#0072C6", "NA" = "grey")
all_drug_list <- unique(all_raw_data$DRUG)
ordered_drugs <- c(drug_list, sort(setdiff(all_drug_list, drug_list)))

All_data_plot <- lapply(ordered_drugs, function(drug) {
  if (drug %in% c("CFZ", "BDQ", "BDQ-M2", "GSK-656", "TBAJ-587", "TBAJ-587-M3", "TBAJ-876", "TBAJ-876-M3")) {
    b <- seq(0, 120, 24)
  } else if (drug == "AMK") {
    b <- seq(0, 24, 2)
  } else if (drug == "TBI-223") {
    b <- seq(0, 24, 4)
  } else {
    b <- seq(0, 72, 6)
  }
  if (drug %in% c("CFZ", "LZD", "DLM", "BDQ", "BDQ-M2", "GSK-286", "SZD", "SZD-M1", "GTX", "TBAJ-876", "TBAJ-876-M3")) {
    yb <- seq(-4, 4, 2)
  } else {
    yb <- seq(-4, 4)
  }
  all_raw_data %>%
    filter(DRUG == drug, DV > 0) %>%
    mutate(LESIONGROUP = ifelse(LESIONGROUP == "Cellular", "Cellular Lesion", LESIONGROUP)) %>%
    filter(LESIONGROUP %in% c("Plasma", "Lung", "Cellular Lesion", "Caseum")) %>%
    mutate(LESIONGROUP = factor(LESIONGROUP, levels = c("Plasma", "Lung", "Cellular Lesion", "Caseum"))) %>%
    ggplot(mapping = aes(x = TAD, y = DV, color = SPECIES, group = SPECIES)) +
    geom_point(size = 1.1, alpha = 0.3) +
    stat_summary_bin(geom = "line", fun = "median", size = 0.3) +
    facet_grid(~LESIONGROUP) +
    ylab(drug) +
    scale_y_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) ifelse(10^x %in% (10^yb), 10^x, NA)),
      labels = trans_format("log10", math_format(10^.x)),
      sec.axis = dup_axis(name = NULL)
    ) +
    scale_x_continuous(breaks = b, limits = c(0, NA)) +
    scale_color_manual(values = SPECIES_color_map) +
    theme_bw() +
    theme(
      legend.background = element_rect("transparent"),
      strip.background = element_rect("transparent"),
      legend.text = element_text(size = 5),
      legend.justification = c(1, 1),
      legend.spacing.y = unit(0, "cm"),
      legend.key.size = unit(0.1, "in"),
      legend.margin = margin(0, 0, 0, 0.01),
      legend.direction = "vertical",
      legend.key.height = unit(0.005, "in"),
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6),
      strip.text.x = element_blank(),
      plot.title = element_text(size = 6, margin = margin(0.1, 0, 0.05, 0, "cm")),
      plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
    )
})

names(All_data_plot) <- ordered_drugs

tiff(file.path(paths$results, "Fig2.tiff"), width = 6.6, height = 8.1, units = "in", res = 300, compression = "lzw")
grid.arrange(
  arrangeGrob(
    grobs = All_data_plot[1:14], ncol = 1,
    top = textGrob("       Plasma           Normal Lung     Cellular Lesion        Caseum", gp = gpar(fontsize = 6), just = 0.52),
    bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 8)),
    left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 8), rot = 90)
  ),
  arrangeGrob(
    grobs = All_data_plot[15:28], ncol = 1,
    top = textGrob("       Plasma           Normal Lung     Cellular Lesion        Caseum", gp = gpar(fontsize = 6), just = 0.52),
    bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 8)),
    left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 8), rot = 90)
  ),
  ncol = 2
)
dev.off()
ppt_add_slide(
  fun = function() {
    grid.arrange(
      arrangeGrob(
        grobs = All_data_plot[1:14], ncol = 1,
        top = textGrob("       Plasma           Normal Lung     Cellular Lesion        Caseum", gp = gpar(fontsize = 6), just = 0.52),
        bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 8)),
        left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 8), rot = 90)
      ),
      arrangeGrob(
        grobs = All_data_plot[15:28], ncol = 1,
        top = textGrob("       Plasma           Normal Lung     Cellular Lesion        Caseum", gp = gpar(fontsize = 6), just = 0.52),
        bottom = textGrob("Time after dose (h)", gp = gpar(fontsize = 8)),
        left = textGrob("Concentration (mg/L or mg/kg)", gp = gpar(fontsize = 8), rot = 90)
      ),
      ncol = 2
    )
  },
  width = 6.6,
  height = 8.1
)
