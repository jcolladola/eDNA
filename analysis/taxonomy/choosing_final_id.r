# run the code with the data in a folder
# Rscript parse_blast.r <blast_output>
library(tidyverse)
params <- commandArgs(TRUE)

target_folder <-  params[1]

taxonomy_cols <- c("phylum", "class", "order", "family", "genus", "species")
output.summary <- function(tibble) {
  tibble %>%
    separate(consensus, into = taxonomy_cols, sep = "%", fill = "right", remove = F) %>%
    mutate(resolution = case_when(
      species  != "NA" ~ "species",
      genus != "NA" ~ "genus",
      family != "NA" ~ "family",
      TRUE ~ "Worse"
    )) |> select(-all_of(taxonomy_cols))
}
ids_outputs <- list.files(path = target_folder, pattern = "ids_threshold_.*\\.tsv", full.names = TRUE) %>%
    set_names(basename(.)) |>
    map_dfr(~read_tsv(.x) |> output.summary(), .id = "file") 





# summary of the resolution levels
ids_outputs %>%
  mutate(threshold = str_extract(file, "(?<=threshold_)\\d+")) %>%
  count(threshold, resolution) %>%
  pivot_wider(names_from = resolution, values_from = n, values_fill = 0)


ids_outputs |>
  mutate(threshold = str_extract(file, "(?<=threshold_)\\d+")) %>%
  pivot_wider(
    id_cols = qseqid,
    names_from = threshold,
    values_from = c(consensus, resolution),
    names_glue = "{.value}_{threshold}"
  ) -> ids_wide

  final_ids <- ids_wide %>%
  mutate(
    final_id = case_when(
      resolution_99 == "species" ~ consensus_99,
      resolution_97 == "species" ~ consensus_97,
      resolution_95 == "species" ~ str_replace(consensus_95, "([^%]*)$", "NA"), # remove species
      resolution_99 == "genus" ~ consensus_99,
      resolution_97 == "genus" ~ consensus_97,
      resolution_95 == "genus" ~ consensus_95,
      resolution_99 == "family" ~ consensus_99,
      resolution_97 == "family" ~ consensus_97,
      resolution_95 == "family" ~ consensus_95,
      TRUE ~ str_replace(consensus_90, "([^%]*)%([^%]*)%([^%]*)$", "NA%NA%NA")

    ),
    reason = case_when(
      resolution_99 == "species" ~ "1",
      resolution_97 == "species" ~ "2",
      resolution_95 == "genus" ~ "3",
      resolution_99 == "genus" ~ "4",
      resolution_97 == "genus" ~ "5",
      resolution_95 == "genus" ~ "6",
      resolution_99 == "family" ~ "7",
      resolution_97 == "family" ~ "8",
      resolution_95 == "family" ~ "9",
      TRUE ~ "10"  # default/fallback
    )
  )

final_ids |> 
  select(qseqid, final_id, reason) |> 
  write_tsv(file.path(target_folder, "final_ids.tsv"))
  
