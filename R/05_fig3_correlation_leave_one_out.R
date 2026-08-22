# Execution order fragile
# Requires: parameters1, lesion_pk_params_rab, lesion_pk_params_hum_updated, drug_list

parms_rab <- parameters1 %>%
  mutate(PC_rab_caseum = ifelse(lesion == "Caseum", PC_rab, NA)) %>%
  arrange(drug) %>%
  group_by(drug) %>%
  fill(PC_rab_caseum, .direction = "downup") %>%
  ungroup() %>%
  drop_na(PC_rab) %>%
  mutate(drug = ifelse(drug == "QBD", "QBS", drug)) %>%
  select(drug, lesion, PC_rab, PC_rab_caseum)

parms_rab$drug <- fct_reorder(parms_rab$drug, parms_rab$PC_rab_caseum, min)
parms_rab$lesion <- factor(parms_rab$lesion, levels = c("Caseum", "Caseous Lesion", "Cavity Wall", "Cellular Lesion", "Normal Lung"))

part_plt <- ggplot(parms_rab, mapping = aes(x = drug, y = PC_rab, shape = lesion, fill = lesion)) +
  geom_point(size = 1.8) +
  geom_point(parms_rab[parms_rab$lesion == "Caseum", ], mapping = aes(x = drug, y = PC_rab), size = 1.8) +
  ylab("Partition coefficient (rabbit)") +
  xlab("") +
  scale_shape_manual(values = c(21, 22, 23, 24, 25)) +
  scale_fill_brewer(palette = "RdBu") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgrey") +
  scale_y_log10() +
  theme_bw(base_size = 10) +
  theme(
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.background = element_rect("transparent"),
    strip.background = element_rect("transparent"),
    legend.text = element_text(size = 10),
    legend.justification = c(0, 1),
    legend.position = c(0.01, 0.99),
    legend.margin = margin(0, 0, 0, 0.1),
    legend.direction = "vertical",
    legend.key.height = unit(0.005, "in"),
    legend.key.size = unit(0.1, "in"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm")
  )

tiff(file.path(paths$results, "Fig3C.tiff"), width = 6.5, height = 3.7, units = "in", res = 300, compression = "lzw")
print(part_plt)
dev.off()
ppt_add_slide(x = part_plt, width = 6.5, height = 3.7)

df_interspecies_corr <- parameters1 %>% filter(SET == 0) %>% drop_na(PC_hu) %>% drop_na(PC_rab) %>% distinct()
DRUG_color_map <- pal_jama("default")(7)
names(DRUG_color_map) <- drug_list

species_corr <- ggscatter(
  df_interspecies_corr, x = "PC_rab", y = "PC_hu",
  add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", cor.coef.size = 3,
  label.x.npc = "right", label.y.npc = "bottom", point = FALSE
) +
  stat_regline_equation(vjust = 0.5, size = 3) +
  geom_abline(slope = 1, linetype = 1) +
  geom_point(df_interspecies_corr, mapping = aes(x = PC_rab, y = PC_hu, shape = lesion, fill = drug), color = "black", size = 1.5) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "") +
  scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) ifelse(10^x %% 1 == 0, 10^x, NA)),
    labels = trans_format("log10", math_format(10^.x)),
    limits = c(NA, 10^2)
  ) +
  scale_x_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) ifelse(10^x %% 1 == 0, 10^x, NA)),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_fill_manual(name = "", values = DRUG_color_map) +
  xlab("Rabbit Partition Coefficient") +
  ylab("Human Partition Coefficient") +
  ggtitle("Interspecies relationship") +
  guides(fill = guide_legend(override.aes = list(shape = 21), size = "none")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    legend.key.height = unit(0.11, "in"),
    legend.key.size = unit(0.1, "in"),
    legend.background = element_rect("transparent"),
    strip.background = element_rect("transparent"),
    strip.text.x = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10, margin = margin(0.1, 0, 0.05, 0, "cm")),
    plot.title = element_text(size = 10, margin = margin(0.1, 0, 0.05, 0, "cm")),
    plot.margin = unit(c(0.1, 0.2, 0.1, 0.15), "cm")
  )

tiff(file.path(paths$results, "Fig3A.tiff"), width = 3, height = 3, units = "in", res = 300, compression = "lzw")
print(species_corr)
dev.off()
ggsave(file.path(paths$results, "Fig3A_species_corr.png"), species_corr, width = 3, height = 3, dpi = 300)
ppt_add_slide(x = species_corr, width = 3, height = 3)

lm_sum <- summary(lm(log10PC_hu ~ log10PC_rab, data = df_interspecies_corr))
final_intercept <- lm_sum[["coefficients"]][1]
final_slope <- lm_sum[["coefficients"]][2]

drug_lst <- unique(df_interspecies_corr$drug)
df <- data.frame()
lst <- list()
for (i in seq_along(drug_lst)) {
  drug_selected <- drug_lst[i]
  remove_one_df <- df_interspecies_corr[df_interspecies_corr$drug != drug_selected, ]
  cortest <- cor.test(remove_one_df$log10PC_hu, remove_one_df$log10PC_rab, method = "pearson")
  lm_sum_loo <- summary(lm(remove_one_df$log10PC_hu ~ remove_one_df$log10PC_rab))

  save_row <- data.frame(drug_missing = drug_selected)
  save_row$intercept <- lm_sum_loo[["coefficients"]][1]
  save_row$slope <- lm_sum_loo[["coefficients"]][2]
  save_row$pearson.cor <- cortest[["estimate"]][["cor"]]
  save_row$p.value <- cortest[["p.value"]]
  df <- rbind(df, save_row)

  plt <- ggscatter(
    remove_one_df, x = "log10PC_rab", y = "log10PC_hu",
    add = "reg.line", conf.int = TRUE, cor.coef = FALSE, cor.method = "pearson",
    label.x.npc = "right", label.y.npc = "bottom", size = 0.5
  ) +
    stat_cor(size = 2) +
    stat_regline_equation(vjust = 0.4, size = 2) +
    geom_point(remove_one_df, mapping = aes(x = log10(PC_rab), y = log10(PC_hu), shape = lesion, fill = drug), size = 2) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "") +
    scale_fill_manual(name = "", values = DRUG_color_map) +
    xlab("Rabbit Partition Coefficient") +
    ylab("Human Partition Coefficient") +
    xlim(-0.75, 2.5) +
    ylim(-0.75, 2.5) +
    ggtitle(paste0(drug_selected, " left out")) +
    guides(fill = guide_legend(override.aes = list(shape = 21), size = "none")) +
    theme_bw(base_size = 6) +
    theme(panel.grid = element_blank(), axis.text = element_text(size = 7), legend.position = "none")
  lst[[length(lst) + 1]] <- plt
}

tiff(file.path(paths$results, "FigS3.tiff"), width = 6.5, height = 4, units = "in", res = 300, compression = "lzw")
grid.arrange(grobs = lst, nrow = 2)
dev.off()
ppt_add_slide(fun = function() grid.arrange(grobs = lst, nrow = 2), width = 6.5, height = 4)

write_csv(df, file.path(paths$results, "FigS3_leave_one_out.csv"))

subset <- df_interspecies_corr %>%
  select(drug, lesion, PC_hu, PC_rab, log10PC_hu, log10PC_rab) %>%
  left_join(df, by = c("drug" = "drug_missing")) %>%
  mutate(
    calc_log10PC_hu_loocv = intercept + (slope * log10PC_rab),
    calc_PC_hu_loocv = round(10^calc_log10PC_hu_loocv, 3),
    calc_log10PC_hu = final_intercept + (final_slope * log10PC_rab),
    calc_PC_hu = round(10^calc_log10PC_hu, 3)
  )

test_set_res <- ggplot() +
  geom_abline(slope = 1, linetype = 1) +
  geom_abline(slope = 1, intercept = 0.5, linetype = 2) +
  geom_abline(slope = 1, intercept = -0.5, linetype = 2) +
  stat_regline_equation(vjust = 0) +
  geom_point(subset, mapping = aes(x = calc_log10PC_hu_loocv, y = log10PC_hu, shape = lesion, fill = drug), color = "black", size = 1.5) +
  scale_shape_manual(values = c(21, 22, 23, 24, 25), name = "") +
  scale_fill_manual(name = "", values = DRUG_color_map) +
  scale_y_continuous(breaks = c(0, 1, 2), labels = math_format(10^.x), limits = c(-1, 2)) +
  scale_x_continuous(breaks = c(0, 1, 2), labels = math_format(10^.x), limits = c(-0.5, 1.5)) +
  xlab("Predicted Partition Coefficient") +
  ylab("Observed Partition Coefficient") +
  ggtitle("Test set results of LOOCV") +
  guides(fill = guide_legend(override.aes = list(shape = 21), size = "none")) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size = 10), legend.position = "none")

tiff(file.path(paths$results, "Fig3B.tiff"), width = 4.5, height = 3, units = "in", res = 300, compression = "lzw")
print(test_set_res)
dev.off()
ppt_add_slide(x = test_set_res, width = 4.5, height = 3)

fig3_top <- ggarrange(species_corr + labs(tags = "A"), test_set_res + labs(tags = "B"), common.legend = TRUE, legend = "right")
fig3_combined <- ggarrange(fig3_top, part_plt + ggtitle("") + labs(tags = "C"), nrow = 2, heights = c(1, 1.5))
tiff(file.path(paths$results, "Fig3.tiff"), width = 6.8, height = 6.8, units = "in", res = 300, compression = "lzw")
print(fig3_combined)
dev.off()
ppt_add_slide(x = fig3_combined, width = 6.8, height = 6.8)

write_csv(subset, file.path(paths$source_data, "predicted_human_partition_coefficients.csv"))
