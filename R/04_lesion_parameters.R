lesion_pk_params_hum_updated <- tribble(
  ~drug, ~hu_KPL_lung, ~hu_KPL_cell, ~hu_KPL_cavity_wall, ~hu_KPL_caseous, ~hu_KPL_caseum, ~hu_PC_lung, ~hu_PC_cell, ~hu_PC_cavity_wall, ~hu_PC_caseous, ~hu_PC_caseum,
  "RIF", 10, 10, 1.98, 0.367, 10, 0.737, 0.4063, 0.614, 0.733, 0.232,
  "INH", 0.146, 0.19, 0.173, 0.168, 0.202, 0.586, 0.522, 0.175, 0.516, 0.719,
  "PZA", 1.83, 10, 1.52, 1.28, 2.12, 0.635, 0.613, 0.72, 0.611, 0.571,
  "MXF", 0.363, 0.77, 0.476, 0.348, 0.354, 4.14, 1.54, 2.33, 2.41, 1.51,
  "CFZ", 4.32, 2.86, 3.21, 2.75, 2.74, 23.27, 14.43, 11.3, 14.3, 2.31,
  "KAN", 0.365, 0.888, 0.186, 0.212, 0.299, 0.555, 0.785, 0.468, 0.503, 0.588,
  "LZD", 0.409, 0.243, 0.305, 0.275, 0.252, 0.755, 0.659, 0.497, 0.177, 0.833
)

lesion_pk_params_rab <- read.table(text = "drug	rab_KPL_lung	rab_KPL_cell	rab_KPL_cavity_wall	rab_KPL_caseous	rab_KPL_caseum	rab_PC_lung	rab_PC_cell	rab_PC_cavity_wall	rab_PC_caseous	rab_PC_caseum
RIF	0.155	10	NA	NA	0.036	2.706	2.907	NA	NA	3.897
INH	0.29	0.06	NA	NA	0.47	1.32	1.15	NA	NA	1.09
PZA	0.544	0.88	1.01	NA	0.616	1.1	0.859	1.49	NA	0.852
MXF	1	0.997	0.726	1.29	10	4.58	8.01	6.14	7.18	3.39
CFZ	10	10	10	10	10	130.4	235.6	NA	NA	2.07
KAN	0.493	0.419	NA	0.448	0.395	0.338	0.454	NA	0.476	0.497
LZD	10	10	10	NA	10	1.04	0.977	1.05	NA	1.05
BDQ	10	0.141	NA	NA	10	133	81.3	NA	NA	4.95
BDQ-M2	0.108	0.034	NA	NA	0.0053	1067	834	NA	NA	31.2
TBAJ-587	10.37	13.16	NA	NA	0.01558	89.76	129.1	NA	NA	7.576
TBAJ-587-M3	10	10	NA	NA	0.6451	973	1172	NA	NA	25.27
TBAJ-876	0.15	0.08	NA	NA	0.03	106.7	145.8	NA	NA	4.82
TBAJ-876-M3	0.04	0.04	NA	NA	0.002	950.9	867	NA	NA	55.47
PMD	2.44	10	10	10	10	4.59	4.56	2.66	4.46	2.16
DLM	2.63	10	NA	0.54	10	8.11	9.18	NA	17.6	1.14
SZD	0.1317	0.2375	NA	NA	0.5334	0.4438	0.5136	NA	NA	0.6626
SZD-M1	0.299	0.228	NA	NA	0.247	0.904	1.01	NA	NA	0.987
TBI-223	10	10	NA	NA	1.72	1.01	1.12	NA	NA	1.49
OPC-167832	10	10	NA	NA	0.393	1.75	1.94	NA	NA	1.60
GSK-656	10	10	NA	NA	10	1.10	1.88	NA	NA	0.77
RPT	0.53	0.82	NA	NA	0.13	1.33	1.99	NA	NA	0.536
GTX	0.699	1.01	NA	NA	10	7.49	8.97	NA	NA	4.62
AMK	0.716	0.385	NA	0.490	0.496	0.437	0.462	NA	0.618	0.927
EMB	8.1	7.52	NA	2.44	2.44	3.90	3.33	NA	2.41	2.41
BTZ-043	0.620	10.8	NA	NA	10	0.831	1.12	NA	NA	0.296
GSK-286	0.39	0.3	NA	NA	0.11	1.62	1.31	NA	NA	1.99
mCLB-073	11.2	2.05	NA	NA	0.19	1.86	1.55	NA	NA	2.08
RBT	0.154	0.209	NA	10	10	29	47.3	NA	51.6	8.94
RZL	10	10	NA	10	10	9.38	23.1	NA	14.4	4.99", header = TRUE)

rab_params_long <- lesion_pk_params_rab %>%
  pivot_longer(-drug, names_to = "term_tag", values_to = "PC_rab") %>%
  filter(grepl("PC", term_tag)) %>%
  mutate(LESIONNAME = case_when(
    grepl("lung", term_tag) ~ "Normal Lung",
    grepl("cell", term_tag) ~ "Cellular Lesion",
    grepl("cavity", term_tag) ~ "Cavity Wall",
    grepl("caseous", term_tag) ~ "Caseous Lesion",
    grepl("caseum", term_tag) ~ "Caseum"
  ))

hum_params_long <- lesion_pk_params_hum_updated %>%
  pivot_longer(-drug, names_to = "term_tag", values_to = "PC_hu") %>%
  filter(grepl("PC", term_tag)) %>%
  mutate(LESIONNAME = case_when(
    grepl("lung", term_tag) ~ "Normal Lung",
    grepl("cell", term_tag) ~ "Cellular Lesion",
    grepl("cavity", term_tag) ~ "Cavity Wall",
    grepl("caseous", term_tag) ~ "Caseous Lesion",
    grepl("caseum", term_tag) ~ "Caseum"
  ))

parameters1 <- full_join(
  hum_params_long %>% select(-term_tag),
  rab_params_long %>% select(-term_tag),
  by = c("drug", "LESIONNAME")
) %>%
  mutate(
    PC_hu = as.numeric(PC_hu),
    PC_rab = as.numeric(PC_rab)
  ) %>%
  mutate(
    log10PC_hu = log10(PC_hu),
    log10PC_rab = log10(PC_rab),
    SET = ifelse(drug %in% drug_list, 0, 1),
    lesion = LESIONNAME,
    drug = ifelse(drug == "GSK-656", "GFB", ifelse(drug == "OPC-167832", "QBS", ifelse(drug == "mCLB-073", "TBD-11", drug)))
  )

write_csv(parameters1, file.path(paths$source_data, "validation_set_partition_coefficients.csv"))
