find_first_existing <- function(candidates) {
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("None of the candidate files exist:\n", paste(candidates, collapse = "\n"))
  }
  hit
}

required_columns <- function(df, cols, df_name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

drug_potency_path <- find_first_existing(c(
  here("data", "derived_datasets", "drug_potency_table.csv"),
  here("drug_potency_table.csv")
))
human_obs_path <- find_first_existing(c(
  here("data", "derived_datasets", "human_lesion_pk.csv"),
  here("human_lesion_pk.csv")
))
rabbit_obs_path <- find_first_existing(c(
  here("data", "derived_datasets", "rabbit_lesion_pk.csv"),
  here("rabbit_lesion_pk.csv")
))

drug_potency_MRT2_edited <- readr::read_csv(drug_potency_path, show_col_types = FALSE)
all_human_obs <- readr::read_csv(human_obs_path, show_col_types = FALSE) %>%
  mutate(METHOD = ifelse(METHOD == "LCM", "LCM/LCMS", "HG/LCMS"))
all_rabbit_obs <- readr::read_csv(rabbit_obs_path, show_col_types = FALSE)

required_columns(all_human_obs, c("DRUG", "CID", "LESIONGROUP"), "all_human_obs")
required_columns(all_rabbit_obs, c("DRUG", "CID", "LESIONGROUP"), "all_rabbit_obs")
find_first_existing <- function(candidates) {
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) {
    stop("None of the candidate files exist:\n", paste(candidates, collapse = "\n"))
  }
  hit
}

required_columns <- function(df, cols, df_name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(df_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

drug_potency_path <- find_first_existing(c(
  here("data", "derived_datasets", "drug_potency_table.csv"),
  here("drug_potency_table.csv")
))

human_obs_path <- find_first_existing(c(
  here("data", "derived_datasets", "human_lesion_pk.csv"),
  here("human_lesion_pk.csv")
))

rabbit_obs_path <- find_first_existing(c(
  here("data", "derived_datasets", "rabbit_lesion_pk.csv"),
  here("rabbit_lesion_pk.csv")
))

drug_potency_MRT2_edited <- readr::read_csv(drug_potency_path, show_col_types = FALSE)
all_human_obs <- readr::read_csv(human_obs_path, show_col_types = FALSE) %>%
  mutate(METHOD = ifelse(METHOD == "LCM", "LCM/LCMS", "HG/LCMS"))
all_rabbit_obs <- readr::read_csv(rabbit_obs_path, show_col_types = FALSE)

required_columns(all_human_obs, c("DRUG", "CID", "LESIONGROUP"), "all_human_obs")
required_columns(all_rabbit_obs, c("DRUG", "CID", "LESIONGROUP"), "all_rabbit_obs")

if (nrow(all_human_obs) == 0 || nrow(all_rabbit_obs) == 0) {
  stop("Imported observation datasets are empty.")
}
