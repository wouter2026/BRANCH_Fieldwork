# -----------------------------
# Setup
# -----------------------------

library(terra)
library(dplyr)
library(tidyr)

# -----------------------------
# Read and inspect file(s)
# -----------------------------

trees_all <- vect("exports/trees_all.gpkg")
names(trees_all)

# -----------------------------
# Create city code + define ID
# -----------------------------
trees_all$city_code <- toupper(substr(trees_all$city, 1, 3))

# Create group ID
group_id <- paste(trees_all$city, trees_all$sample_round)

# Generate sequence per group
trees_all$tree_num <- ave(
  seq_len(nrow(trees_all)),
  group_id,
  FUN = function(x) sprintf("%05d", seq_along(x))
)

# Final ID
trees_all$ID_sampling <- paste(
  trees_all$city_code,
  trees_all$sample_round,
  trees_all$tree_num,
  sep = "_"
)



# -----------------------------
# clean columns
# -----------------------------

# Define columns to keep
cols_keep <- c("id", "species_raw", "genus", "ID_sampling","sample_round")

# Extract attributes + geometry preserved
trees_cleaned <- trees_all[, cols_keep]

# -----------------------------
# Add placeholder fields
# -----------------------------

trees_cleaned$Date <- as.Date(NA)
trees_cleaned$Name_s <- NA_character_
trees_cleaned$Species <- NA_character_
trees_cleaned$DBH <- NA_real_
trees_cleaned$Crown_Dia_NS <- NA_real_
trees_cleaned$Crown_Dia_EW <- NA_real_
trees_cleaned$Height <- NA_real_
trees_cleaned$Height_base_crown <- NA_real_
trees_cleaned$LAI <- NA_real_
trees_cleaned$n_veg_layer <- NA_integer_

trees_cleaned$Context_under <- NA_character_
trees_cleaned$Context_around <- NA_character_
trees_cleaned$Context_pit <- NA_character_

trees_cleaned$Comments <- NA_character_
trees_cleaned$Photo <- NA_character_

# Reorder

trees_cleaned <- trees_cleaned[, c(
  "id", "species_raw", "genus", "ID_sampling", "sample_round",
  "Date", "Name_s", "Species", "DBH",
  "Crown_Dia_NS", "Crown_Dia_EW",
  "Height", "Height_base_crown",
  "LAI", "n_veg_layer",
  "Context_under", "Context_around", "Context_pit",
  "Comments", "Photo"
)]

trees_all$city_code <- toupper(substr(trees_all$city, 1, 3))

# Create group ID
group_id <- paste(trees_all$city, trees_all$sample_round)

# Generate sequence per group
trees_all$tree_num <- ave(
  seq_len(nrow(trees_all)),
  group_id,
  FUN = function(x) sprintf("%05d", seq_along(x))
)

# Final ID
trees_all$ID_sampling <- paste(
  trees_all$city_code,
  trees_all$sample_round,
  trees_all$tree_num,
  sep = "_"
)

# -----------------------------
# Split by sample_round
# -----------------------------

rounds <- unique(trees_cleaned$sample_round)

trees_split <- lapply(rounds, function(r) {
  trees_cleaned[trees_cleaned$sample_round == r, ]
})

names(trees_split) <- paste0("trees_round_", rounds)

# -----------------------------
# export files
# -----------------------------

# Create function
export_gpkg <- function(vect, name, folder = "exports") {
  if (!dir.exists(folder)) dir.create(folder)
  
  file <- file.path(folder, paste0(name, ".gpkg"))
  
  # Write as GeoPackage
  writeVector(vect, file, overwrite = TRUE)
  
  message("Exported spatial file: ", file)
}

# Export  datasets with all trees
export_gpkg(trees_all, "trees_all")
export_gpkg(trees_cleaned, "trees_all_cleaned")

# export split datasets

for (name in names(trees_split)) {
  export_gpkg(trees_split[[name]], name)
}